library(tidyverse)
library(xlsx)
bwa.comm <- '~/software/bwa/bwa mem -M -t 10 data/genome/BWA/GCF_000470375.1_Vir5453'
genomeAnnotationFile <- 'data/genome/GCF_000470375.1_Vir5453_genomic.gbff'
sampleDataDir <- 'data/sampleSeqData'
outputDir <- 'ESSENTIALS_input'
serverURL <- 'http://bushmanlab.org/data/export/everett/StGeme/'


sampleNames <- 
  list('GTSP3171' = 'Challenge_1.1', 'GTSP3172' = 'Challenge_1.2',
       'GTSP3173' = 'Challenge_2.1', 'GTSP3174' = 'Challenge_2.2',
       'GTSP3175' = 'Challenge_3.1', 'GTSP3176' = 'Challenge_3.2',
       'GTSP3177' = 'Control_1.1',   'GTSP3178' = 'Control_1.2',
       'GTSP3179' = 'Control_2.1',   'GTSP3180' = 'Control_2.2',
       'GTSP3181' = 'Control_3.1',   'GTSP3182' = 'Control_3.2')

d <- tibble(file = list.files(path = sampleDataDir, pattern = 'GTSP'),
            GTSP = str_extract(file, 'GTSP\\d+'),
            direction = ifelse(grepl('virusReads', file), 'transposon', 'breakPoint'))

d$sampleName <- unname(sapply(d$GTSP, function(x) sampleNames[[x]]))

config <- bind_rows(lapply(split(d, d$sampleName), function(x){
  if(n_distinct(x$direction) != 2) stop('Direction error.')
  comm <- paste0(bwa.comm, ' ', sampleDataDir, '/', x[which(x$direction == 'breakPoint'),]$file, ' ',
                 sampleDataDir, '/', x[which(x$direction == 'transposon'),]$file, ' > ESSENTIALS_input/', x[1,]$sampleName, '.sam')
  system(comm)
  tibble(link = paste0(serverURL, x[1,]$sampleName, '.sam'),
         barcode = 'N',
         transposon_end = 'N',
         sampletype = ifelse(grepl('Control', x[1,]$sampleName, ignore.case = TRUE), 'control', 'target'),
         library = 'lib1',
         sampleformat = 'sam',
         compression = 'none')
}))

write.table(config, file = 'ESSENTIALS_input/config.tsv', col.names = TRUE, row.names = FALSE, sep = '\t')
write.xlsx(config, file = 'ESSENTIALS_input/config.xlsx', col.names = TRUE, row.names = FALSE)
system(paste0('cp ', genomeAnnotationFile, ' ESSENTIALS_input/genBankAnnoations'))
system('scp ESSENTIALS_input/* microb120:/media/lorax/data/export/everett/StGeme/')
