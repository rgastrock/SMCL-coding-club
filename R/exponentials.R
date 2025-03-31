


# function optimization -----

exponentialModel <- function(par, timepoints) {
  
  if (length(timepoints) == 1) {
    timepoints <- c(0:(timepoints-1))
  }
  
  # print(par)
  # print(timepoints)
  
  output = par['N0'] - ( par['N0'] * (1-par['lambda'])^timepoints )
  
  # print(output)
  
  return(data.frame(trial=timepoints,
                    output=output))
  
}


exponentialMSE <- function(par, signal, timepoints=c(0:(length(signal)-1)) ) {
  
  MSE <- mean((exponentialModel(par, timepoints)$output - signal)^2, na.rm=TRUE)
  
  return( MSE )
  
}


exponentialFit <- function(participants=NULL, df, timepoints=length(signal), gridpoints=4, gridfits=4, asymptoteRange=NULL) {
  
  
  # select relevant participants
  subdf <- NA
  
  for (participant in participants) {
    
    pdf <- df[which(df$participant == participant),]
    
    if (is.data.frame(subdf)) {
      subdf <- rbind(subdf, pdf)
    } else {
      subdf <- pdf
    }
    
  }
  
  agdf <- aggregate(reachdeviation_deg ~ trial_num, data=subdf, FUN=mean, na.rm=TRUE)
  signal <- agdf$reachdeviation_deg
  timepoints=length(signal)
  
  
  
  # set the search grid:
  parvals <- seq(1/gridpoints/2,1-(1/gridpoints/2),1/gridpoints)
  
  if (is.null(asymptoteRange)) {
    # set a wiiiiide range... especially for single participants, the range may or may not work depending on how noisy their data is
    asymptoteRange <- c(-1,2)*max(abs(signal), na.rm=TRUE)
  }
  
  
  searchgrid <- expand.grid('lambda' = parvals,
                            'N0'     = parvals * diff(asymptoteRange) + asymptoteRange[1] )
  
  # evaluate starting positions:
  MSE <- apply(searchgrid, FUN=exponentialMSE, MARGIN=c(1), signal=signal, timepoints=timepoints)
  
  
  
  
  lo <- c(0,asymptoteRange[1])
  hi <- c(1,asymptoteRange[2])
  
  
  
  
  
  
  # run optimx on the best starting positions:
  allfits <- do.call("rbind",
                     apply( data.frame(searchgrid[order(MSE)[1:gridfits],]),
                            MARGIN=c(1),
                            FUN=optimx::optimx,
                            fn=exponentialMSE,
                            method     = 'L-BFGS-B',
                            lower      = lo,
                            upper      = hi,
                            timepoints = timepoints,
                            signal     = signal
                     ) )
  
  # pick the best fit:
  win <- allfits[order(allfits$value)[1],]
  
  winpar <- unlist(win[1:2])
  
  # return the best parameters:
  return(winpar)
  
}


# fit exponential to learning curves -----


groupLearningExponentials <- function() {
  
  # loop through groups
  for (group in c('control', 'cursorjump', 'handview')) {
    
    # load the group reach training data:
    df <- read.csv(sprintf('data/%s/%s_training_reachdevs.csv', group, group),
                   stringsAsFactors = FALSE)
    
    # fit exponential:
    exp_par <- exponentialFit( participants = unique(df$participant), df=df)
    
    # print best parameters for now:
    print(exp_par)
    
  }
  
}


bootstrapExponentials <- function(bootstraps=200) {
  
  # set up a cluster:
  ncores <- parallel::detectCores()
  clust  <- parallel::makeCluster(max(c(1,floor(ncores*0.80))))
  
  parallel::clusterEvalQ(cl=clust, source('R/exponentials.R'))
  
  # loop through groups
  for (group in c('control', 'cursorjump', 'handview')) {
    
    cat('working on group:', group, '\n')
    
    # load the group reach training data:
    df <- read.csv(sprintf('data/%s/%s_training_reachdevs.csv', group, group),
                   stringsAsFactors = FALSE)
    
    # create matrix of randomly sampled participants to bootstrap across participants:
    participants <- unique(df$participant)
    BSparticipants <- matrix( sample(participants,
                                     size=bootstraps*length(participants),
                                     replace=TRUE),
                              nrow=bootstraps)
    
    # fit an exponential to bootstrapped sets of participants' data:
    a <- parallel::parApply(cl = clust,
                            X = BSparticipants,
                            MARGIN = 1,
                            FUN = exponentialFit,
                            df = df)
    
    # write the fits to a file:
    outdf <- as.data.frame(t(a))
    write.csv(outdf,
              file=sprintf('data/%s/%s_expfits.csv',group,group),
              row.names = FALSE)
    
  }
  
  # stop the cluster and free the cores for other tasks:
  parallel::stopCluster(clust)
  
}


# analyse bootstrapped exponentials -----

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
