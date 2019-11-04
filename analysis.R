library(gt23)
library(parallel)
library(GenomicRanges)
library(tidyverse)
library(ggrepel)
library(ggbeeswarm)
source('analysis.lib.R')
load('data/sites.RData')

# select(subset(sites, start >= 1900944 - 5 & start <= 1900944 + 5), start, posid, reads, uniqueSample, estAbund) %>% arrange(uniqueSample, start) %>% select(-start)



cluster <- makeCluster(30)
clusterExport(cluster, 'sites')
pos <- subset(sites, strand == '+')
pos$s <- ntile(1:nrow(pos), 30)

duelDetections <- 
  bind_rows(parLapply(cluster, split(pos, pos$s), function(a){
       library(dplyr)
       s <- split(a, paste(a$uniqueSample, a$posid))
       bind_rows(lapply(s, function(x){
              m <- subset(sites, start >= x$start-5 & start <= x$start & uniqueSample == x$uniqueSample & strand == '-')
              if(nrow(m) == 1){
                return(tibble(pos = paste0(x$uniqueSample, '~', x$posid), neg =  paste0(x$uniqueSample, '~', m$posid), maxAbund = max(x$estAbund, m$estAbund)))
              } else if (nrow(m) > 1){
                return(tibble(pos = paste0(x$uniqueSample, '~', x$posid), neg =  'multiple', maxAbund = NA))
              } else {
                return(tibble(pos = paste0(x$uniqueSample, '~', x$posid), neg = 'none', maxAbund = NA))
              }
       }))
}))

duelDetections <- subset(duelDetections, ! is.na(maxAbund))

sites$s <- paste0(sites$uniqueSample, '~', sites$posid)
clusterExport(cluster, 'sites')
duelDetections$s <- ntile(1:nrow(duelDetections), 30)

mergedSites <- bind_rows(parLapply(cluster, split(duelDetections, duelDetections$s), function(z){
                  bind_rows(lapply(1:nrow(z), function(x){
                    x <- z[x,]
                    a <- subset(sites, s == x$pos)
                    b <- subset(sites, s == x$neg)
                    a$mergedSiteIDs <- list(c(a$s, b$s))
                    a$reads <- sum(a$reads, b$reads)
                    a$estAbund <- x$maxAbund
                    a$strand <- '*'
                    a$start <- floor(mean(a$start, b$start))
                    a$end <- a$start
                    a
                  }))
              }))

mergedSitesIDs <- unlist(mergedSites$mergedSiteIDs)
mergedSites$mergedSiteIDs <- NULL
sites <- bind_rows(mergedSites, subset(sites, ! s %in% mergedSitesIDs))



report <- list()
sites$patient   <- sites$subject
sites$GTSP      <- sites$sample
sites$refGenome <- 'Kingella_kingae_KKKWG1'
sites$sample2 <- paste0(sites$subject, '.', sites$replicate)

sites <- makeGRangesFromDataFrame(sites, keep.extra.columns = TRUE) %>%
         annotateIntSites() %>%
         data.frame() %>%
         group_by(sample2) %>%
         mutate(nSampleSites = n_distinct(posid),
                nSampleCells = sum(estAbund),
                siteCounts = 1 / nSampleSites[1]) %>%
         ungroup()


sites.tbl <- group_by(sites, sample2) %>%
             summarise(insertions = n_distinct(posid)) %>%
             ungroup()



        
m1 <- bind_rows(lapply(split(sites, sites$sample2), function(x){
        tibble(a = seq(1, 2000000, by = 10000),
               b = seq(10000, 2000000, by = 10000),
               sample = x$sample2[1]) %>%
        mutate(s =  1:n()) %>%
        rowwise() %>%
        mutate(n = sum(subset(x, start >= a & start < b)$siteCounts)) %>%
        ungroup()
      }))
     


m2 <- bind_rows(lapply(split(sites, sites$sample2), function(x){
        tibble(a = seq(1, 2000000, by = 10000),
               b = seq(10000, 2000000, by = 10000),
               sample = x$sample2[1]) %>%
        mutate(s =  1:n()) %>%
        rowwise() %>%
        mutate(n = sum(subset(x, start >= a & start < b)$estAbund) / x$nSampleCells[1]) %>%
        ungroup()
     }))


sampleLevels <- rev(c('Control_1.1', 'Control_1.2', 'Challenge_1.1', 'Challenge_1.2',
                      'Control_2.1', 'Control_2.2', 'Challenge_2.1', 'Challenge_2.2',
                      'Control_3.1', 'Control_3.2', 'Challenge_3.1', 'Challenge_3.2'))

m1$sample <- factor(as.character(m1$sample), levels = sampleLevels)
m2$sample <- factor(as.character(m2$sample), levels = sampleLevels)


genomeInsertionPlot <- 
  ggplot(m1, aes(s, sample, fill = n)) +
  geom_tile() +
  scale_fill_gradient(name = 'Normalized insertion counts', low = 'white', high = "black",  
                      na.value = "grey50") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) +
  theme(text = element_text(size=14)) +
  labs(x = '10KB genomic blocks', y = 'Experiments') +
  theme(legend.position="bottom") +
  guides(fill=guide_colourbar(barwidth=25,label.position="top"))


genomeAbundancePlot <- 
  ggplot(m2, aes(s, sample, fill = n)) +
  geom_tile() +
  scale_fill_gradient(name = 'Normalized cell counts', low = 'white', high = "black",  
                      na.value = "grey50") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) +
  theme(text = element_text(size=14)) +
  labs(x = '10KB genomic blocks', y = 'Experiments') +
  theme(legend.position="bottom") +
  guides(fill=guide_colourbar(barwidth=25,label.position="top"))


# Here we normalize the number of sites within each gene by 
# dividing the number of insertions in gene X by the total
# number of insertions in the sample.
d1 <- expandData(group_by(sites, sample2) %>%
      mutate(nSampleSites = n_distinct(posid)) %>%
      ungroup() %>%
      filter(inFeature == TRUE) %>%
      group_by(sample2, nearestFeature) %>%
      summarise(nSitesNorm = n_distinct(posid) / nSampleSites[1]) %>%
      ungroup()) 

d1.result <- collapseRepsGeneEnrichment(d1)


d2 <- expandData(group_by(sites, sample2) %>%
                   mutate(nSampleCells = sum(estAbund)) %>%
                   ungroup() %>%
                   filter(inFeature == TRUE) %>%
                   group_by(sample2, nearestFeature) %>%
                   summarise(nSitesNorm = sum(estAbund) / nSampleCells[1]) %>%
                   ungroup())  

d2.result <- collapseRepsGeneEnrichment(d2)


save.image(file = 'data/analysis.RData')


library(xlsx)
write.xlsx2(d1.result$enrichmentTable, 'siteCountApproachTable.xlsx', col.names = TRUE, row.names = FALSE)
write.xlsx2(d2.result$enrichmentTable, 'abundanceApproachTable.xlsx', col.names = TRUE, row.names = FALSE)


