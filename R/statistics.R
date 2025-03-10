getANOVAdata <- function() {
  
  aov_data <- NA
  
  for (group in c('control','cursorjump','handview')) {
    
    gdf <- read.csv(sprintf('data/%s/%s_training_reachdevs.csv', group, group),
                    stringsAsFactors = FALSE)
    
    gdf$block <- NA
    
    gdf$block[which(gdf$trial_num %in% c(1,2,3))] <- 1
    gdf$block[which(gdf$trial_num %in% c(4,5,6))] <- 2
    gdf$block[which(gdf$trial_num %in% c(76:90))] <- 3
    
    bgdf <- aggregate(reachdeviation_deg ~ block + participant,
                      data = gdf,
                      FUN = mean,
                      na.rm=TRUE)
    
    bgdf$group <- group
    
    if (is.data.frame(aov_data)) {
      aov_data <- rbind(aov_data, bgdf)
    } else {
      aov_data <- bgdf
    }
    
  }
  
  aov_data <- aov_data[which(aov_data$participant != '6468c3'),]
  
  aov_data$group <- as.factor(aov_data$group)
  aov_data$block <- as.factor(aov_data$block)
  aov_data$participant <- as.factor(aov_data$participant)
  
  
  return(aov_data)
  
}

learningCurveANOVA <- function() {
  
  aov_data <- getANOVAdata()
  
  aovm <- afex::aov_ez( id = 'participant',
                        dv = 'reachdeviation_deg',
                        data = aov_data,
                        between = c('group'),
                        within = c('block'))
  
  print(aovm)
  
}


learningCurveBayes <- function() {
  
  aov_data <- getANOVAdata()
  
  # baov <- BayesFactor::anovaBF(reachdeviation_deg ~ group + block,
  #                              data = aov_data,
  #                              whichRandom = "participant",
  #                              progress=FALSE)
  # print(baov)
  
  
  block3 <- aov_data[which(aov_data$block == 3),]
  
  baov <- BayesFactor::anovaBF(reachdeviation_deg ~ group,
                               data = block3,
                               whichRandom = "participant",
                               progress=FALSE)
  
  print(baov)
  
}