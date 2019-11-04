# Conversion table for interal IDs to human readable IDs.
sampleNames <-  list('GTSP3171' = 'Challenge_1.1', 'GTSP3172' = 'Challenge_1.2',
                     'GTSP3173' = 'Challenge_2.1', 'GTSP3174' = 'Challenge_2.2',
                     'GTSP3175' = 'Challenge_3.1', 'GTSP3176' = 'Challenge_3.2',
                     'GTSP3177' = 'Control_1.1',   'GTSP3178' = 'Control_1.2',
                     'GTSP3179' = 'Control_2.1',   'GTSP3180' = 'Control_2.2',
                     'GTSP3181' = 'Control_3.1',   'GTSP3182' = 'Control_3.2')

expandData <- function(d){
  # Here we cast out the data into a wide format in order to store zeros for instances 
  # where there were no insertions in one or more replicates after which we melt the data.
  reshape2::dcast(select(d, sample2, nearestFeature, nSitesNorm), sample2~nearestFeature, value.var = 'nSitesNorm', fill = 0) %>%
  reshape2::melt(id.vars = 'sample2', variable.name = 'nearestFeature', value.name = 'nSitesNorm')
}


collapseRepsGeneEnrichment <- function(d){

  d.table <- spread(d, key = sample2, value = nSitesNorm)

  d$sample <- sub('\\.\\d$', '', d$sample2)
  
  d.wide <- reshape2::dcast(select(d, sample2, nearestFeature, nSitesNorm), sample2~nearestFeature, value.var = 'nSitesNorm')
  rownames(d.wide) <- d.wide$sample2
  d.wide$sample2 <- NULL
  pca <- prcomp(d.wide, scale = TRUE, center = TRUE)
  plotData  <- data.frame(s = row.names(pca$x), x = pca$x[,1], y = pca$x[,2], z = pca$x[,3])
  plotData$sampleType <- ifelse(grepl('Control', plotData$s), 'Control', 'Challenge')
  
  
  set.seed(46)
  PCA_plot <-
    ggplot(plotData, aes(x=x, y=y, color = sampleType, label = s)) +
    theme_bw() +
    geom_point(size = 2, stroke = 1.5) +
    geom_text_repel(point.padding = 0.2) +
    scale_color_manual(values = c('gray30', 'gray70')) +
    labs(x = paste0('PC1 (', sprintf("%.2f", summary(pca)$importance[3,][1] * 100), '%)'),
         y = paste0('PC2 (', sprintf("%.2f", (summary(pca)$importance[3,][2] - summary(pca)$importance[3,][1])  * 100), '%)'))
  
  
  d.averagedReps <- group_by(d, sample, nearestFeature) %>%
                    summarise(mean_nSitesNorm = mean(nSitesNorm), nReplicates = n()) %>%
                    ungroup()
  
  
  
  d.wide <- reshape2::dcast(select(d.averagedReps, sample, nearestFeature, mean_nSitesNorm), sample~nearestFeature, value.var = 'mean_nSitesNorm')
  rownames(d.wide) <- d.wide$sample
  d.wide$sample <- NULL
  pca <- prcomp(d.wide, scale = TRUE, center = TRUE)
  plotData  <- data.frame(s = row.names(pca$x), x = pca$x[,1], y = pca$x[,2], z = pca$x[,3])
  plotData$sampleType <- ifelse(grepl('Control', plotData$s), 'Control', 'Challenge')
  
  set.seed(46)
  PCA_plot_reps_averaged <-
    ggplot(plotData, aes(x=x, y=y, color = sampleType, label = s)) +
    theme_bw() +
    geom_point(size = 2, stroke = 1.5) +
    geom_text_repel(point.padding = 0.2) +
    scale_color_manual(values = c('gray30', 'gray70')) +
    labs(x = paste0('PC1 (', sprintf("%.2f", summary(pca)$importance[3,][1] * 100), '%)'),
         y = paste0('PC2 (', sprintf("%.2f", (summary(pca)$importance[3,][2] - summary(pca)$importance[3,][1])  * 100), '%)'))
  
  d.replicateDifference <- bind_rows(lapply(split(d, paste(d$sample, d$nearestFeature)), function(x){
    tibble(sample = x[1,]$sample, diff = x[1,]$nSitesNorm - x[2,]$nSitesNorm)
  }))

  d.differencePlot <- 
    ggplot(mapping = aes(d.replicateDifference$sample, d.replicateDifference$diff)) + 
    geom_quasirandom(alpha=.2) +
    theme_bw() +
    labs(x = 'Experiment', y = 'Normalized number of insertions') +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_flip()

  
  getGeneDesc <- function(x){
    desc <- paste0(sapply(unlist(strsplit(as.character(x), ',')), 
                    function(g) gt23::Kingella_kingae_KKKWG1.refSeqGenesGRanges[match(g, gt23::Kingella_kingae_KKKWG1.refSeqGenesGRanges$name2),]$name), collapse = ', ')
    gsub('\\n\\s*', ' ', desc)
  }
  
  
  d.geneEnrcihments <- bind_rows(lapply(split(d, as.character(d$nearestFeature)), function(x){
    control   <- x[which(grepl('Control', x$sample)),]$nSitesNorm
    challenge <- x[-which(grepl('Control', x$sample)),]$nSitesNorm
    
    result <- tryCatch({
      p <- t.test(control, challenge)$p.value
    },  warning = function(war) {
      return(p)
    }, error = function(err) {
      return(NA)
    })
    
    tibble(nearestFeature = x$nearestFeature[1], 
           pVal = result,
           higherInChallenge = ifelse(mean(challenge) > mean(control), TRUE, FALSE),
           geneDesc = getGeneDesc(nearestFeature))
  })) %>%  arrange(pVal) %>% mutate(pVal.adj = p.adjust(pVal))
  
  
  
  d.geneEnrcihments.averagedReps <- bind_rows(lapply(split(d.averagedReps, as.character(d.averagedReps$nearestFeature)), function(x){
    control   <- x[which(grepl('Control', x$sample)),]$mean_nSitesNorm
    challenge <- x[-which(grepl('Control', x$sample)),]$mean_nSitesNorm
  
    result <- tryCatch({
      p <- t.test(control, challenge)$p.value
    },  warning = function(war) {
      return(p)
    }, error = function(err) {
      return(NA)
    })
  
    tibble(nearestFeature = x$nearestFeature[1], 
           pVal_avgReps = result,
           higherInChallenge_avgReps = ifelse(mean(challenge) > mean(control), TRUE, FALSE))
  })) %>% arrange(pVal_avgReps) %>% mutate(pVal.adj_avgReps = p.adjust(pVal_avgReps))


  d.table <- left_join(d.table, d.geneEnrcihments, by = 'nearestFeature') 
  d.table <- left_join(d.table, d.geneEnrcihments.averagedReps, by = 'nearestFeature') 
  
  return(list(diffPlot = d.differencePlot, enrichmentTable = d.table, PCA = PCA_plot, PCA_avg = PCA_plot_reps_averaged))
}


heatmap_dims <- function(p) {
  .x <- as.character(p$mapping$x)
  .y <- as.character(p$mapping$y)
  
  .x <- .x[-grep('~', .x)]
  .y <- .y[-grep('~', .y)]
  
  ncols <- length(unique(p$data[[.x]]))
  nrows <- length(unique(p$data[[.y]]))
  return(list(ncols=ncols, nrows=nrows))
}

make_square <- function(p, fudge=1) {
  dims <- heatmap_dims(p)
  p + ggplot2::theme(aspect.ratio = (dims$nrows/dims$ncols)*fudge)
}









