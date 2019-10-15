library(yaml)
library(ShortRead)
library(tidyverse)
library(parallel)
library(GenomicRanges)
library(gintools)

config  <- read_yaml('~/releases/TransposonMapping_StGeme/StGeme.config')

source(file.path(config$visaHome, 'visa.lib.R'))
config$startTime <- proc.time()

samples <- read_delim(config$sampleConfigFile, delim  = ',', col_names = TRUE, col_types = cols())

# Create unique sample identifiers. 
if(any(grepl('~|\\|', paste(samples$subject, samples$sample, samples$replicate)))) 
  stop('Error -- tildas (~) are reserved characters and can not be used in the subject, sample, or replicate sample configuration columns.')


# Create unique sample identifiers (subject ~ sample id ~ replicate).
samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)


# Reverse compliment the provided index1 barcodes if requested.
if(config$indexReads.rc) samples$index1Seq <- as.character(reverseComplement(DNAStringSet(samples$index1Seq)))


# Chunk the read data.
if(! dir.exists(config$outputDir)) dir.create(config$outputDir)
dir.create(file.path(config$outputDir, 'tmp'))
dir.create(file.path(config$outputDir, 'readIDs'))
dir.create(file.path(config$outputDir, 'readsRemoved'))
dir.create(file.path(config$outputDir, 'seqChunks'))
dir.create(file.path(config$outputDir, 'sampleReads'))

cluster <- makeCluster(config$demultiplexing.CPUs)
clusterExport(cluster, c('config', 'samples'))

invisible(parLapply(cluster, 
                    list(c(config$index1ReadsFile, config$sequence.chunk.size, 'index1Reads', file.path(config$outputDir, 'seqChunks')),
                         c(config$breakReadsFile,  config$sequence.chunk.size, 'breakReads',  file.path(config$outputDir, 'seqChunks')),
                         c(config$virusReadsFile,  config$sequence.chunk.size, 'virusReads',  file.path(config$outputDir, 'seqChunks'))), 
                 function(x){
                   library(ShortRead)
                   source(file.path(config$visaHome, 'visa.lib.R'))
                   createSeqChunks(x[[1]], x[[2]], x[[3]], x[[4]])
                }))


if(config$correctGolayIndexReads){
  invisible(parLapply(cluster, list.files(file.path(config$outputDir, 'seqChunks'), full.names = TRUE, pattern = 'index'), function(x){
  #invisible(lapply(list.files(file.path(config$outputDir, 'seqChunks'), full.names = TRUE, pattern = 'index'), function(x){
    library(ShortRead)
    library(dplyr)
    source(file.path(config$visaHome, 'visa.lib.R'))
    golayCorrection(x)
  }))
}




write(date(), file =  file.path(config$outputDir, 'log'), append = FALSE)
invisible(parLapply(cluster, list.files(file.path(config$outputDir, 'seqChunks'), pattern = 'virusReads', full.names = TRUE), function(f){
#invisible(lapply(list.files(file.path(config$outputDir, 'seqChunks'), pattern = 'virusReads', full.names = TRUE), function(f){
  library(ShortRead)
  library(tidyverse)
  source(file.path(config$visaHome, 'visa.lib.R'))
  
  # Capture the chunk identifier.
  chunk.n <- unlist(str_match_all(f, '\\.(\\d+)$'))[2]
  
  
  # Create a local instance of the config object to store the chunk start time.
  te <- paste(sprintf("%.1f", (proc.time() - config$startTime)[3] / 60), 'minutes elapsed since start.')
  config$chunkStartTime <- proc.time() 
  write(paste0('Starting sequence chunk ', chunk.n, '. ', te), file = file.path(config$outputDir, 'log'), append = TRUE)
  
  chunkMsg <- function(chunk.n, config, msg){
    te <- paste(sprintf("%.1f", (proc.time() - config$chunkStartTime)[3] / 60), 'minutes elapsed since sequence chunk start.')
    write(paste0('Sequence chunk ', chunk.n, '.\t', msg, '\t', te), file = file.path(config$outputDir, 'log'), append = TRUE)
  }
  
  
  # Read in sequence chunks.
  virusReads  <- readFastq(f)
  breakReads  <- readFastq(sub('virusReads', 'breakReads',  f))
  index1Reads <- readFastq(sub('virusReads', 'index1Reads', f))
  
  
  # Trim reads by base call quality scores using a sliding window approach.
  preTrimReadIDs <-  sub('\\s+.+$', '', as.character(virusReads@id))
  if(length(preTrimReadIDs) > 0) writeLines(preTrimReadIDs, file.path(config$outputDir, 'readIDs', paste0(chunk.n, '.readIDs')))
  
  virusReads <- trimTailw(virusReads, 2, config$sequence.trim.qualCode, 5)
  breakReads <- trimTailw(breakReads, 2, config$sequence.trim.qualCode, 5)
  
  if(length(virusReads) == 0 | length(breakReads) == 0){
    chunkMsg(chunk.n, config, 'No reads remaining after quality trimming.')
    return()
  }

  
  # Convert reads to DNAstring sets in order to conserve memory and more easily access read sequences.
  index1Reads <- shortRead2DNAstringSet(index1Reads)
  virusReads  <- shortRead2DNAstringSet(virusReads)
  breakReads  <- shortRead2DNAstringSet(breakReads)
  
  
  # Sync becuase reads may of been lost in the previous base call quality filter.
  reads <- syncReads(index1Reads, virusReads, breakReads)
  index1Reads <- reads[[1]];  virusReads  <- reads[[2]];  breakReads  <- reads[[3]]
  
  
  # Report reads removed from quality trimming.
  ids <- preTrimReadIDs[! preTrimReadIDs %in% names(virusReads)]
  if(length(ids) > 0) writeLines(ids, file.path(config$outputDir, 'readsRemoved', paste0('qualTrim.', chunk.n)))
  
  if(length(virusReads) == 0){
    chunkMsg(chunk.n, config, 'No reads remaining after quality trimming and read synchronization.')
    return()
  }
  
  # Loop through samples in sample data file.
  invisible(lapply(1:nrow(samples), function(r){
    r <- samples[r,]
    
    # Trim over-read sequences.
    virusReads <- trimOverReadSeq(virusReads, r$virusRead.overReadSeq)
    breakReads <- trimOverReadSeq(breakReads, r$breakRead.overReadSeq)
    
    
    # Ensure that reads are long enough for upcoming sample specific tests.
    preReadSizeCheck1 <- names(virusReads)
    virusReads <- virusReads[width(virusReads) >= (nchar(r$virusLTRSeq) + config$trimmedRead.minLength)]
    breakReads <- breakReads[width(breakReads) >= (nchar(r$breakReadLinkerSeq) + config$trimmedRead.minLength)]
    
    
    # Sync becuase reads may of been lost in the previous length filter.
    reads <- syncReads(index1Reads, virusReads, breakReads)
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]];  breakReads  <- reads[[3]]
    
    
    # Report reads removed from min read size filter1.
    ids <- preReadSizeCheck1[! preReadSizeCheck1 %in% names(virusReads)]
    if(length(ids) > 0) writeLines(ids, file.path(config$outputDir, 'readsRemoved', paste0('readSizeCheck1.', r$uniqueSample, '.', chunk.n)))
    
    if(length(virusReads) == 0){
      chunkMsg(chunk.n, config, paste0('(', r$uniqueSample, ') No reads remaining after first read size filter.'))
      return()
    }
    
    # test
    # M03249:18:000000000-CN6GB:1:2108:16355:3459 ???
    # browser()
    
    # Create barcode demultiplexing vector.
    v1 <- rep(TRUE, length(virusReads))
    if('index1Reads.maxMismatch' %in% names(config)){
      v1 <- vcountPattern(r$index1Seq, index1Reads, max.mismatch = config$index1Reads.maxMismatch) > 0
    }
    
    
    # Create break read linker barcode demultiplexing vector.
    v2 <- rep(TRUE, length(breakReads))
    if('breakReads.linkerBarcode.maxMismatch' %in% names(config)){
      testSeq <- substr(r$breakReadLinkerSeq, r$breakReadLinkerBarcode.start, r$breakReadLinkerBarcode.end)
      v2 <- vcountPattern(testSeq, subseq(breakReads, r$breakReadLinkerBarcode.start, r$breakReadLinkerBarcode.end), max.mismatch = config$breakReads.linkerBarcode.maxMismatch) > 0
    }
    
    
    # Test to see if any reads demultiplex to this row of the sample table and then subset reads to this sample.
    i <- base::intersect(which(v1), which(v2))
    if(length(i) == 0){
      chunkMsg(chunk.n, config, paste0('No reads demultiplexed for sample: ', r$uniqueSample, '.'))
      return()
    }
    
    reads <- syncReads(index1Reads[i], virusReads[i], breakReads[i])
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
    if(length(virusReads) == 0){
      chunkMsg(chunk.n, config, paste0('(', r$uniqueSample, ') No reads remaining after post demultiplexing read synchronization.'))
      return()
    }
    

    # Test vector read alignemnts to the vector sequences.
    vectorReadIDs <- getVectorReadIDs(virusReads, config, r$vectorSeqFile)
    
    if(length(vectorReadIDs) > 0){
      virusReadsVector <- virusReads[names(virusReads) %in% vectorReadIDs]
      breakReadsVector <- breakReads[names(breakReads) %in% vectorReadIDs]
    
      names(virusReadsVector) <- paste0(names(virusReadsVector), '|', r$uniqueSample)
      names(breakReadsVector) <- paste0(names(breakReadsVector), '|', r$uniqueSample)
      
      writeFasta(virusReadsVector, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.virusReadsVector.', chunk.n, '.fasta')))
      writeFasta(breakReadsVector, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.breakReadsVector.', chunk.n, '.fasta')))
    }
    
    # Filter virus reads against read ids returned by vectorReadIDs().
    virusReads <- virusReads[! names(virusReads) %in% vectorReadIDs]
    
    
    # Sync reads. 
    reads <- syncReads(index1Reads, virusReads, breakReads)
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
    if(length(virusReads) == 0){
      chunkMsg(chunk.n, config, paste0('(', r$uniqueSample, ') No reads remaining after vector read removal and read synchronization.'))
      return()
    }

    
    # Test the start of virus reads.
    v1 <- rep(TRUE, length(virusReads))
    if('virusReads.startTest.maxMismatch' %in% names(config)){
      testSeq <- substr(r$virusLTRSeq, 1,  config$virusReads.startTest.length)
      v1 <- vcountPattern(testSeq, subseq(virusReads, 1, config$virusReads.startTest.length), max.mismatch = config$virusReads.startTest.maxMismatch) > 0
    }
    
    
    # Test the entire LTR sequence.
    v2 <- rep(TRUE, length(virusReads))
    if('virusReads.fullTest.maxMismatch' %in% names(config)){
      v2 <- vcountPattern(r$virusLTRSeq, subseq(virusReads, 1, nchar(r$virusLTRSeq)), max.mismatch = config$virusReads.fullTest.maxMismatch) > 0
    }
    
    
    # Test break read common linker.
    v3 <- rep(TRUE, length(virusReads))
    if('breakReads.linkerCommon.maxMismatch' %in% names(config)){
      testSeq <- substr(r$breakReadLinkerSeq, r$breakReadLinkerCommon.start, r$breakReadLinkerCommon.end)
      v3 <- vcountPattern(testSeq, subseq(breakReads, r$breakReadLinkerCommon.start,r$breakReadLinkerCommon.end), max.mismatch = config$breakReads.linkerCommon.maxMismatch) > 0
    }
    
    
    # Test to see which reads pass the last three filters.
    i <- base::Reduce(base::intersect, list(which(v1), which(v2), which(v3)))
    if(length(i) == 0){
      chunkMsg(chunk.n, config, paste0('(', r$uniqueSample, ') No reads remaining after user specified read tests and read synchronization.'))
      return()
    }
    
  
    # Subset and sync reads.
    reads <- syncReads(index1Reads[i], virusReads[i], breakReads[i])
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
  
      
    # Trim leading adapter sequences.
    virusReads <- trimLeadingSeq(virusReads, r$virusLTRSeq)
    breakReads <- trimLeadingSeq(breakReads, r$breakReadLinkerSeq)
    
    
    # Select reads which have the minimum read lengths post trimming.
    preFilterReads <- names(virusReads)
    virusReads <- virusReads[width(virusReads) >= config$trimmedRead.minLength]
    breakReads <- breakReads[width(breakReads) >= config$trimmedRead.minLength]
    
    
    # Report reads removed from min read size filter2.
    ids <- preFilterReads[! preFilterReads %in% names(virusReads)]
    if(length(ids) > 0) writeLines(ids, file.path(config$outputDir, 'readsRemoved', paste0('readSizeCheck2.', r$uniqueSample, '.', chunk.n)))
    
    
    # Sync reads after selecting for reads with min. length post-trimming.
    reads <- syncReads(index1Reads, virusReads, breakReads)
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
    if(length(virusReads) == 0){
      chunkMsg(chunk.n, config, paste0('(', r$uniqueSample, ') No reads remaining after second read size filter.'))
      return()
    }
    
  
    # Write out final reads. Add sample names to read IDs.
    names(breakReads) <- paste0(names(breakReads), '|', r$uniqueSample)
    names(virusReads) <- paste0(names(virusReads), '|', r$uniqueSample)
    
    writeFasta(breakReads,  file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.breakReads.', chunk.n, '.fasta')))
    writeFasta(virusReads,  file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.virusReads.', chunk.n, '.fasta')))
    
    chunkMsg(chunk.n, config, paste0('(', r$uniqueSample, ') completed with ', length(virusReads), ' reads.'))
  }))
  
  chunkMsg(chunk.n, config, 'Chunk completed.')
}))
write(date(), file = file.path(config$outputDir, 'log'), append = TRUE)
stopCluster(cluster)


# Collate demultiplexed reads using file name snibets.
collateSampleReads('virusReads')
collateSampleReads('breakReads')
collateSampleReads('virusReadsVector')
collateSampleReads('breakReadsVector')

system(paste('cat ', paste0(list.files(file.path(config$outputDir, 'sampleReads'), full.names = TRUE, pattern = 'virusReadsVector'), collapse = ' '), '> ', file.path(config$outputDir, 'vectorVirusReads.fasta')))
system(paste('cat ', paste0(list.files(file.path(config$outputDir, 'sampleReads'), full.names = TRUE, pattern = 'breakReadsVector'), collapse = ' '), '> ', file.path(config$outputDir, 'vectorBreakReads.fasta')))


# Organize read files by reference genome and read source (virusRead or breakRead).
# This will allow samples to be aligned to different genomes defined in the samples table.
d <- tibble(file = list.files(file.path(config$outputDir, 'sampleReads'), pattern = '\\.virusReads\\.|\\.breakReads\\.'),
            uniqueSample = unlist(lapply(strsplit(file, '\\.'), function(x) paste0(x[1:(length(x) - 2)], collapse = '.'))),
            source = ifelse(grepl('virusRead', file), 'virusReads', 'breakReads')) %>%
     left_join(select(samples, uniqueSample, refGenomeBLATdb), by = 'uniqueSample')



# Align the sample FASTA files defined in the previous data frame to their respective genomes. 
# This is done by grouping reads by type (virusReads / breakReads) and reference geneome 
# and using parLapply() to disrubte read chunks across a number of CPUs.

alignments <- bind_rows(lapply(split(d, paste(d$source, d$refGenomeBLATdb)), function(x){
  
  f <- tmpFile()
  system(paste('cat ', paste0(file.path(config$outputDir, 'sampleReads', x$file), collapse = ' '), ' > ', 
               paste0(file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))))
  
  fasta <- readFasta(file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))
  file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))
  
  db <- x$refGenomeBLATdb[1]
  cluster <- makeCluster(config$genomAlignment.CPUs)
  clusterExport(cluster, c('config', 'db'), envir = environment())
  
  r <- bind_rows(parLapply(cluster, split(fasta, ceiling(seq_along(fasta)/config$alignment.chunk.size)), function(y){
         library(ShortRead)
         library(tidyverse)
         source(file.path(config$visaHome, 'visa.lib.R'))
         alignReads.BLAT(y, db)
       }))
  
  stopCluster(cluster)
  
  r$source <- x$source[1]
  r
}))


# Apply alignment filters.
alignments <- 
  filter(alignments, 
         alignmentPercentID   >= config$alignment.genome.minPercentID,
         percentQueryCoverage >= config$alignment.genome.minPercentQueryCoverage,
         qStart <= 5,
         tNumInsert <= 1, 
         qNumInsert <= 1,
         tBaseInsert <= 2,
         qBaseInsert <= 2) 


# M03249:18:000000000-CN6GB:1:1101:25986:21460


# Remove multihits
v <- table(subset(alignments, source == 'virusReads')$qName)
alignments <- subset(alignments, ! qName %in% names(v[v > 1]))


# Rename alignment column headers by appending read sources.
alignments <- lapply(split(alignments, alignments$source), function(x){
    names(x) <- paste0(names(x), '.', x$source[1])
    x
})

# Combine read alignments into rows, extract sample names from read ids and collapse identical alignments. 
frags <- 
  left_join(alignments[["virusReads"]], alignments[["breakReads"]], by = c('qName.virusReads' = 'qName.breakReads')) %>%
  select(qName.virusReads, strand.virusReads, strand.breakReads, tName.virusReads, tName.breakReads,
         tStart.virusReads, tStart.breakReads, tEnd.virusReads, tEnd.breakReads) %>%
  drop_na() %>%
  mutate(uniqueSample = unlist(lapply(str_split(qName.virusReads, '\\|'), '[[', 2 ))) %>%
  group_by(uniqueSample, strand.virusReads, strand.breakReads, tName.virusReads, tName.breakReads,
           tStart.virusReads, tStart.breakReads, tEnd.virusReads, tEnd.breakReads) %>%
  summarise(readCount = n(), readIDs = list(qName.virusReads)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(fragStart = ifelse(strand.virusReads == '+', tStart.virusReads + 1, tStart.breakReads + 1),
         fragEnd   = ifelse(strand.virusReads == '+', tEnd.breakReads + 1,   tEnd.virusReads + 1),
         fragTest  = ifelse(strand.virusReads == '+', tStart.virusReads < tEnd.breakReads, tStart.breakReads < tEnd.virusReads),  
         fragWidth = (fragEnd - fragStart) + 1) %>%
  ungroup() %>%
  filter(fragTest == TRUE,
         strand.breakReads != strand.virusReads & 
         tName.breakReads == tName.virusReads & 
         fragWidth <= config$fragments.maxLength & 
         fragWidth >= config$fragments.minLength) %>%
  mutate(fragID = paste0(uniqueSample, ';', tName.virusReads, ';', strand.virusReads, ';', fragStart, ';', fragEnd)) %>%
  group_by(fragID) %>%
  summarise(reads = sum(readCount)) %>%
  ungroup() %>%
  separate(fragID, c('uniqueSample', 'seqnames', 'strand', 'start', 'end'), sep = ';') %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


frags$subject <- unlist(lapply(strsplit(frags$uniqueSample, '~'), '[[', 1))

cluster <- makeCluster(config$demultiplexing.CPUs)

frags.std <- unlist(GRangesList(parLapply(cluster, split(frags, frags$subject), function(x) gintools::standardize_sites(x, counts.col = 'reads')))) 

frags.std <- 
  unlist(GRangesList(parLapply(cluster, split(frags.std, frags.std$uniqueSample), function(x) gintools::refine_breakpoints(x, counts.col = 'reads')))) %>%
  data.frame() %>%
  group_by(uniqueSample, subject,seqnames, start, end, strand) %>%
  summarise(reads = sum(reads)) %>%
  ungroup()

stopCluster(cluster)
  

sites <- 
  mutate(frags.std, posid = paste0(seqnames, strand, ifelse(strand == '+', start, end))) %>%
  group_by(uniqueSample, posid) %>%
  mutate(reads = sum(reads), estAbund = n()) %>%
  slice(1) %>%
  mutate(start = ifelse(strand == '+', start, end), end = start) %>%
  ungroup() %>%
  group_by(uniqueSample) %>%
  mutate(relAbund = estAbund / sum(estAbund)) %>%
  ungroup()



d <- strsplit(sites$uniqueSample, '~')
sites$uniqueSample <- NULL
sites$subject   <- unlist(lapply(d, '[[', 1))
sites$sample    <- unlist(lapply(d, '[[', 2))
sites$replicate <- unlist(lapply(d, '[[', 3))


save(sites, file = file.path(config$outputDir, 'sites.RData'))
