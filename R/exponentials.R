
#function optimization----

exponentialModel <- function(par, timepoints, mode='learning', setN0=NULL) {
  
  if (length(timepoints) == 1) {
    timepoints <- c(0:(timepoints-1))
  }
  
  if (is.numeric(setN0)) {
    par['N0'] = setN0
  }
  
  if (mode == 'learning') {
    output = par['N0'] - ( par['N0'] * (1-par['lambda'])^timepoints )
  }
  if (mode == 'washout') {
    output = par['N0'] * (par['lambda'])^timepoints
  }
  
  return(data.frame(trial=timepoints,
                    output=output))
  
}


exponentialMSE <- function(par, signal, timepoints=c(0:(length(signal)-1)), mode='learning', setN0=NULL) {
  
  MSE <- mean((Reach::exponentialModel(par, timepoints, mode=mode, setN0=setN0)$output - signal)^2, na.rm=TRUE)
  
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
  MSE <- apply(searchgrid, FUN=Reach::exponentialMSE, MARGIN=c(1), signal=signal, timepoints=timepoints)
  
  
  
  
  lo <- c(0,asymptoteRange[1])
  hi <- c(1,asymptoteRange[2])
  
  
  
  
  
  
  # run optimx on the best starting positions:
  allfits <- do.call("rbind",
                     apply( data.frame(searchgrid[order(MSE)[1:gridfits],]),
                            MARGIN=c(1),
                            FUN=optimx::optimx,
                            fn=Reach::exponentialMSE,
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

#fit exponential learning curves----

groupLearningExponentials <- function(){
  
  for (group in c('control', 'cursorjump', 'handview')){
    
    df <- read.csv(sprintf('data/%s/%s_training_reachdevs.csv', group, group), stringsAsFactors = F)
    
    # adf <- aggregate(reachdeviation_deg ~ trial_num, data=df, FUN=mean, na.rm=T)
    
    exp_par <- exponentialFit(participants = unique(df$participant), df=df)
    print(exp_par)
    
  }
  
  
  
}

bootStrapExponentials <- function(bootstraps=200) {
  
  # set up a cluster:
  ncores <- parallel::detectCores()
  clust  <- parallel::makeCluster(max(c(1,floor(ncores*0.80))))
  # clust  <- parallel::makeCluster(2)
  
  parallel::clusterEvalQ(cl=clust, expr="source('R/exponentials.R')")
  
  
  # loop through groups
  for (group in c('control', 'cursorjump', 'handview')) {
    
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
    write.csv(outdf, file=sprintf('data/%s/%s_expfits.csv',group,group))
    
  }
  
  # stop the cluster and free the cores for other tasks:
  parallel::stopCluster(clust)
  
  
}





















