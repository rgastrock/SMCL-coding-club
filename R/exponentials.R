
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


exponentialFit <- function(signal, timepoints=length(signal), mode='learning', gridpoints=11, gridfits=10, setN0=NULL,asymptoteRange=NULL) {
  
  # set the search grid:
  parvals <- seq(1/gridpoints/2,1-(1/gridpoints/2),1/gridpoints)
  
  if (is.null(asymptoteRange)) {
    # set a wiiiiide range... especially for single participants, the range may or may not work depending on how noisy their data is
    asymptoteRange <- c(-1,2)*max(abs(signal), na.rm=TRUE)
  }
  
  # define the search grid:
  # if (is.numeric(setN0)) {
  #   searchgrid <- expand.grid('lambda' = parvals)
  #   lo <- c(0)
  #   hi <- c(1)
  # }
  if (is.null(setN0)) {
    searchgrid <- expand.grid('lambda' = parvals,
                              'N0'     = parvals * diff(asymptoteRange) + asymptoteRange[1] )
    lo <- c(0,asymptoteRange[1])
    hi <- c(1,asymptoteRange[2])
  } else {
    searchgrid <- expand.grid('lambda' = parvals,
                              'N0'     = setN0)
    lo <- c(0,setN0)
    hi <- c(1,setN0)
  }
  # evaluate starting positions:
  MSE <- apply(searchgrid, FUN=exponentialMSE, MARGIN=c(1), signal=signal, timepoints=timepoints, mode=mode, setN0=setN0)
  
  # if (is.null(setN0)) {
  #   X <- data.frame(searchgrid[order(MSE)[1:gridfits],])
  # } else {
  #   X <- data.frame('lambda'=searchgrid[order(MSE)[1:gridfits],])
  # }
  
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
                            signal     = signal,
                            mode       = mode,
                            setN0      = setN0 ) )
  
  # pick the best fit:
  win <- allfits[order(allfits$value)[1],]
  
  if (is.null(setN0)) {
    winpar <- unlist(win[1:2])
  } else {
    winpar <- c( 'lambda' = unlist(win[1]), 
                 'N0'     = setN0)
    names(winpar) <- c('lambda', 'N0')
  }
  
  # return the best parameters:
  return(winpar)
  
}

#fit exponential learning curves----

learningExponentials <- function(){
  
  for (group in c('control', 'cursorjump', 'handview')[1]){
    
    df <- read.csv(sprintf('data/%s/%s_training_reachdevs.csv', group, group), stringsAsFactors = F)
    
    adf <- aggregate(reachdeviation_deg ~ trial_num, data=df, FUN=mean, na.rm=T)
    
    exp_par <- exponentialFit(signal = adf$reachdeviation_deg)
    
    
  }
  
  
  
}