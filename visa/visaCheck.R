library(ShortRead)

# Why did 'M03249:18:000000000-CN6GB:1:2108:16355:3459' demultiplex to GTSP3174 ?
# linker barcode: CGGCTTACAATTCCTGCGAC
# ITR start: CAACCTGTTA
# bar code:  CGGGATCAAATT (AATTTGATCCCG)                                              
# A,GTSP3174,1,1,20,CGGCTTACAATTCCTGCGACNNNNNNNNNNNNCTCCGCTTAAGGGACT,CGGGATCAAATT,CAACCTGTTA,
#                   CGGCTTACAATTCCTGCGAC                             AATTTGATCCCG
# 
# CAACCTGTTAAGTAATCGTAAACGATGAAAATAATAGGCTGCGGTATTTTTATTAACACCTACTAACTCTGCTGCTGTTCTTGCGGTTACACCTGTGACAAATAGCTCAATGAGTTTGTTTTGTTTATACTGACTTAGACGAC
# ~/software/ncbi-blast-2.7.1+/bin/blastn -db data/genome/BLAST/GCF_000470375.1_Vir5453_genomic.fna  -query  test.ff  -out test.blast
# 
# M03249:18:000000000-CN6GB:1:1101:20903:906
# CAACCTGTTA AATATACTCTGTTCTCGCTGCTTTTTTATTGCTTATTTTTTTACAGCACACTTTTCAAGTCCCTTAAGCGGAGTTAGTTAAGTTGGTCGCAGGAATTGTAAGCCGACAATTACCATAGCGTCAGTCCTGGTGT
# 
# from sampleReads:
# >M03249:18:000000000-CN6GB:1:1101:20903:9069|A~GTSP3174~1
# AATATACTCTGTTCTCGCTGCTTTTTTATTGCTTATTTTTTTACAGCACACTTTTCA
# 
# ~/software/ncbi-blast-2.7.1+/bin/blastn -db data/genome/BLAST/GCF_000470375.1_Vir5453_genomic.fna  -query  test.ff
# (nothing)
# 
# match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts  tStarts
#         match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# 55      1       0       0       1       1       0       0       -       test    57      0       57      chr1    1992631 1717636 1717692 2       21,35,  0,22,   1717636,1717657,
# 55      1       0       0       1       1       0       0       -       test    57      0       57      chr1    1992631 168474  168530  2       21,35,  0,22,   168474,168495,
# 
# 

g <- readFasta('data/genome/GCF_000470375.1_Vir5453_genomic.fna')
as.character(subseq(sread(g), 168474, 168530))
"CTGAAAAGTGTGCTGTAAAAAATAAGCAATAAAAAAGCAGCGAAAACAGAGTATATT"
  TGAAAAGTGTGCTGAAAAAAATAAGCAATAAAAAAGCAGCGAGAACAGAGTATATT

testSample <- 'visaOutput/sampleReads/A~GTSP3174~1.virusReads.fasta'

I1 <- readFastq('/data/sequencing/Illumina-archive/191002_M03249_0018_000000000-CN6GB/Data/Intensities/BaseCalls/Undetermined_S0_L001_I1_001.fastq.gz')
R1 <- readFastq('/data/sequencing/Illumina-archive/191002_M03249_0018_000000000-CN6GB/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz')
R2 <- readFastq('/data/sequencing/Illumina-archive/191002_M03249_0018_000000000-CN6GB/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz')

truncated.R2.ids <- sub('\\s.+$', '', as.character(R2@id))

test <- readFasta(testSample)

testSample.ids   <- unlist(lapply(strsplit(as.character(test@id), '\\|'), '[[', 1))

sample.I1 <- I1[truncated.R2.ids %in% testSample.ids]
sample.R1 <- R1[truncated.R2.ids %in% testSample.ids]
sample.R2 <- R2[truncated.R2.ids %in% testSample.ids]

I1.tab <- data.frame(sort(table(as.character(sread(sample.I1))), decreasing = TRUE))
I1.tab$editDist <- adist(I1.tab$Var1, 'AATTTGATCCCG')
write.table(I1.tab, file = 't', quote = FALSE) 


sort(table(as.character(subseq(sread(sample.R1), 1, 20))), decreasing = TRUE)

