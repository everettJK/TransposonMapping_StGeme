

golayCorrection <- function(x){
  fq <- readFastq(x)
  fq.ids <- sub('\\s+.+$', '', as.character(fq@id))
  
  writeFasta(fq, file = file.path(paste0(x, '.fasta')))
  
  system(paste(config$command.python2, file.path(config$visaHome, 'golayCorrection', 'processGolay.py'), file.path(paste0(x, '.fasta'))))
  invisible(file.remove(c(paste0(x, '.fasta'), x)))
  
  corrected <- readFasta(file.path(paste0(x, '.fasta.corrected')))
  
  a <- left_join(tibble(id = fq.ids, seq = as.character(sread(fq))),
                 tibble(id = as.character(corrected@id), seq2 = as.character(sread(corrected))), by = 'id')
  a <- rowwise(a) %>% mutate(editDist = as.integer(adist(seq, seq2))) %>% ungroup()
  
  i <- which(a$editDist <= 2)
  a[i,]$seq <- a[i,]$seq2
  
  if(! all(fq.ids == a$id)) stop('There was an ordering error during the Golay correction step.')
  
  o <- paste0('@', as.character(fq@id), '\n', a$seq, '\n+\n', as.character(quality(fq@quality)))
  fileConn <- file(x)
  writeLines(o, fileConn)
  close(fileConn)
}



createSeqChunks <- function(f, chunkSize, label, ouputDir){
  strm <- FastqStreamer(f, n = as.integer(chunkSize))
  n <- 1
  repeat {
    fq <- yield(strm)
    if(length(fq) == 0) break
    writeFastq(fq, file = file.path(ouputDir, paste0(label, '.', n)), compress = FALSE)
    n <- n + 1
  }
}


collateSampleReads <- function(label){
  v <- tibble(file = list.files(file.path(config$outputDir, 'tmp'), pattern = paste0('\\.', label, '\\.')),
              sample = unlist(lapply(strsplit(file, paste0('\\.', label, '\\.')), '[[', 1)))
  invisible(lapply(split(v, v$sample), function(x){
    system(paste('cat ', paste0(file.path(config$outputDir, 'tmp', x$file), collapse = ' '), ' > ', 
                 paste0(file.path(config$outputDir, 'sampleReads', paste0(x$sample[1], '.', label, '.fasta')))))
  }))
}



alignReads.BLAT <- function(x, db){
  f <- tmpFile()
  writeFasta(x, file = file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))
  
  comm <- paste0(config$command.blat, ' ', db, ' ', 
                 file.path(config$outputDir, 'tmp', paste0(f, '.fasta')), ' ', 
                 file.path(config$outputDir, 'tmp', paste0(f, '.psl')),  
                 ' -tileSize=11 -stepSize=9 -minIdentity=90 -out=psl -noHead')
  system(comm)
  file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))
  
  if(file.exists(file.path(config$outputDir, 'tmp', paste0(f, '.psl')))){
    b <- parseBLAToutput(file.path(config$outputDir, 'tmp', paste0(f, '.psl')))
    file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.psl')))
    b
  } else {
    return(tibble())
  }
}

# BLAT output parser.
parseBLAToutput <- function(f){
  if(! file.exists(f) | file.info(f)$size == 0) return(tibble())
  b <- read_delim(f, delim = '\t', col_names = FALSE, col_types = cols())
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  
  b$queryPercentID       <- (b$matches/b$qSize)*100
  b$tAlignmentWidth      <- (b$tEnd - b$tStart) + 1
  b$queryWidth           <- (b$qEnd - b$qStart) + 1
  b$alignmentPercentID   <- (b$matches/b$tAlignmentWidth)*100
  b$percentQueryCoverage <- (b$queryWidth/b$qSize)*100
  b$qStarts              <- as.character(b$qStarts)
  b$tStarts              <- as.character(b$tStarts)
  b
}


tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }


shortRead2DNAstringSet <- function(x){
  r <- x@sread
  names(r) <- sub('\\s+.+$', '', as.character(x@id))
  r
}


waitForFile <- function(f, seconds = 1){
  repeat
  {
    if(file.exists(f)) break
    Sys.sleep(seconds)
  }
  
  return(TRUE)
}


syncReads <-function(...){
  arguments <- list(...)
  n <- Reduce(base::intersect, lapply(arguments, names))
  lapply(arguments, function(x){ x <- x[names(x) %in% n]; x[order(names(x))]})
}


syncShortReads <- function(x){
  n <- Reduce(base::intersect, list(as.character(ShortRead::id(x$I1)), as.character(ShortRead::id(x$R1)), as.character(ShortRead::id(x$R2))))
  x$I1 <- x$I1[as.character(ShortRead::id(x$I1)) %in% n]
  x$R1 <- x$R1[as.character(ShortRead::id(x$R1)) %in% n]
  x$R2 <- x$R2[as.character(ShortRead::id(x$R2)) %in% n]
  x$R1 <- x$R1[match(as.character(ShortRead::id(x$R1)), as.character(ShortRead::id(x$I1)))]
  x$R2 <- x$R2[match(as.character(ShortRead::id(x$R2)), as.character(ShortRead::id(x$I1)))]
  message('Reads remaining: ', length(x$I1))
  x
}



trimLeadingSeq <- function(x, seq){
  f <- tmpFile()
  
  writeFasta(x,  file = file.path(config$outputDir, 'tmp', f))
  
  system(paste0(config$command.cutadapt3, ' -f fasta  -e 0.15 -g ', seq, ' --overlap 2 ', 
                file.path(config$outputDir, 'tmp', f), ' > ', 
                file.path(config$outputDir, 'tmp', paste0(f, '.out'))),
         ignore.stderr = TRUE)
  
  s <- readFasta(file.path(config$outputDir, 'tmp', paste0(f, '.out')))
  x <- s@sread
  names(x) <- as.character(s@id)
  file.remove(file.path(config$outputDir, 'tmp', f))
  file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.out')))
  x
}


trimOverReadSeq <- function(x, seq){
  f <- tmpFile()
  
  writeFasta(x,  file = file.path(config$outputDir, 'tmp', f))
  
  system(paste0(config$command.cutadapt3, ' -f fasta  -e 0.15 -a ', seq, ' --overlap 2 ', 
                file.path(config$outputDir, 'tmp', f), ' > ', 
                file.path(config$outputDir, 'tmp', paste0(f, '.out'))),
         ignore.stderr = TRUE)
  
  s <- readFasta(file.path(config$outputDir, 'tmp', paste0(f, '.out')))
  x <- s@sread
  names(x) <- as.character(s@id)
  file.remove(file.path(config$outputDir, 'tmp', f))
  file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.out')))
  x
}






getVectorReadIDs <- function(reads, config, vectorFile){
  d <- group_by(tibble(ids = names(reads), 
                       read = as.character(reads)), read) %>%
    summarise(id.list = paste0(ids, collapse=',')) %>%
    ungroup() %>%
    mutate(id = paste0('s', 1:n()),
           qlength = nchar(read))
  
    f <- tmpFile()
    write(paste0('>', d$id, '\n', d$read), file = file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))

    system(paste0(config$command.blastBin, '/makeblastdb -in ', vectorFile, ' -dbtype nucl -out ', file.path(config$outputDir, 'tmp', f)), ignore.stderr = TRUE)
   
    waitForFile(file.path(file.path(config$outputDir, 'tmp', paste0(f, '.nin'))))
    
    system(paste0(config$command.blastBin, '/blastn -word_size 10 -evalue 10 -outfmt 6 -query ',
                  file.path(config$outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
                  file.path(config$outputDir, 'tmp', f),
                  ' -num_threads 5 -out ', file.path(config$outputDir, 'tmp', paste0(f, '.blast'))), 
           ignore.stdout = TRUE, ignore.stderr = TRUE)
  
    waitForFile(file.path(config$outputDir, 'tmp', paste0(f, '.blast')), seconds = 1)
  
    if(file.info(file.path(config$outputDir, 'tmp', paste0(f, '.blast')))$size == 0) return(character(length = 0))
    
    b <- read.table(file.path(config$outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    b <- left_join(b, d, by = c('qname' = 'id'))
    b$pcoverage <- ((b$qend - b$qstart) / b$qlength)*100
  
    b <- subset(b, pident >= config$alignment.vector.minPercentID & pcoverage >= config$alignment.vector.minPercentQueryCoverage & gapopen <= 1)

    invisible(file.remove(list.files(file.path(config$outputDir, 'tmp'), pattern = f, full.names = TRUE)))
    unlist(strsplit(b$id.list, ','))
}

