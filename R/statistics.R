

# learning curve statistics ----

## regular NHST / Anova ----

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
  
  # aov_data <- aov_data[which(aov_data$participant != '6468c3'),]
  
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

## Bayesian statistics ---- 

learningCurveBayes <- function() {
  
  aov_data <- getANOVAdata()
  
  # baov <- BayesFactor::anovaBF(reachdeviation_deg ~ group * block + participant,
  #                              data = aov_data,
  #                              whichRandom = "participant",
  #                              progress=FALSE)
  # print(baov)
  
  
  block3 <- aov_data[which(aov_data$block == 3),]
  
  baov <- BayesFactor::anovaBF(reachdeviation_deg ~ group,
                               data = block3,
                               # whichRandom = "participant",
                               progress=FALSE)
  
  print(baov)
  
}



## analyse bootstrapped exponentials -----

analyseExponentials <- function() {
  
  # loop through groups
  for (group in c('control', 'cursorjump', 'handview')) {
    
    # load the bootstrapped exponential fits:
    expfits <- read.csv(sprintf('data/%s/%s_expfits.csv', group, group),
                        stringsAsFactors = FALSE)
    
    # calculate the 95% CI for lambda and N0
    lambdaCI <- quantile(expfits$lambda, c(0.025, 0.50, 0.975))
    N0CI     <- quantile(expfits$N0, c(0.025, 0.50, 0.975))
    
    cat(sprintf('lambda 95%% CI for %s: mean=%0.3f, range: %0.3f, %0.3f\n', group, lambdaCI[2], lambdaCI[1], lambdaCI[3]))
    
  }
  
}


# no-cursor statistics ----

## regular NHST / Anova ----

doNoCursorANOVA <- function() {
  
  df <- getNoCursorData()
  
  aovm <- afex::aov_ez( id = 'participant',
                        dv = 'reachdeviation_deg',
                        data = df,
                        between = c('group'),
                        within = c('strategy'))
  
  print(aovm)
  
  cellmeans <- emmeans::emmeans(aovm, specs=c('group','strategy'))
  
  contrasts <- list('control'                = c(1, 0, 0,-1, 0, 0),
                    'NS_control_cursorjump'  = c(1,-1, 0, 0, 0, 0),
                    'NS_control_handview'    = c(1, 0,-1, 0, 0, 0),
                    'WS_control_cursorjump'  = c(0, 0, 0, 1,-1, 0),
                    'WS_control_handview'    = c(0, 0, 0, 1, 0,-1),
                    'NS_cursorjump_handview' = c(0, 1,-1, 0, 0, 0))
  
  emmeans::contrast(cellmeans, contrasts, adjust='sidak')
  
}
