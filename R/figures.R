

# generic plotting functions -----

setupFigureFile <- function(target='inline',width=8,height=6,dpi=300,filename) {
  
  if (target == 'pdf') {
    pdf(file   = filename, 
        width  = width, 
        height = height)
  }
  if (target == 'svg') {
    svglite::svglite( filename = filename,
                      width = width,
                      height = height,
                      fix_text_size = FALSE) 
    # fix_text_size messes up figures on my machine... 
    # maybe it's better on yours?
  }
  if (target == 'png') {
    png( filename = filename,
         width = width*dpi,
         height = height*dpi,
         res = dpi
    )
  }
  if (target == 'tiff') {
    tiff( filename = filename,
          compression = 'lzw',
          width = width*dpi,
          height = height*dpi,
          res = dpi
    )
  }
}



getColors <- function() {
  
  cols.op <- c(rgb(255, 147, 41,  255, max = 255), # orange:  21, 255, 148
               rgb(229, 22,  54,  255, max = 255), # red:    248, 210, 126
               rgb(207, 0,   216, 255, max = 255), # pink:   211, 255, 108
               rgb(127, 0,   216, 255, max = 255), # violet: 195, 255, 108
               rgb(0,   19,  136, 255, max = 255)) # blue:   164, 255, 68
  
  cols.tr <- c(rgb(255, 147, 41,  32,  max = 255), # orange:  21, 255, 148
               rgb(229, 22,  54,  32,  max = 255), # red:    248, 210, 126
               rgb(207, 0,   216, 32,  max = 255), # pink:   211, 255, 108
               rgb(127, 0,   216, 32,  max = 255), # violet: 195, 255, 108
               rgb(0,   19,  136, 32,  max = 255)) # blue:   164, 255, 68
  
  cols <- list()
  cols$op <- cols.op
  cols$tr <- cols.tr
  
  return(cols)
  
}

# training plots -----

plotTraining <- function(target='inline') {
  
  
  setupFigureFile(target=target,
                  width = 8,
                  height=6,
                  dpi=300,
                  sprintf('doc/fig2_training.%s', target))
  
  groups <- c('control', 'cursorjump', 'handview')
  
  colors <- getColors()
  
  plot(NA, NA,
       main='learning', xlab='trial', ylab='reach deviation [°]',
       xlim=c(0,91), ylim=c(-10,40),
       ax=F, bty='n')
  
  lines(x=c(0,91),
        y=c(30,30),
        col='#999999',
        lw=2,
        lty=2)
  lines(x=c(0,91),
        y=c(0,0),
        col='#999999',
        lw=2,
        lty=2)
  
  for (group_no in c(1:length(groups))) {
    
    group <- groups[group_no]
    
    
    df <- read.csv(sprintf('data/%s/%s_training_reachdevs.csv', group, group),
                   stringsAsFactors = F)
    
    
    CI <- aggregate(reachdeviation_deg ~ trial_num,
                    data = df,
                    FUN = Reach::getConfidenceInterval)
    
    lo <- CI$reachdeviation_deg[,1]
    hi <- CI$reachdeviation_deg[,2]
    
    polygon( x = c( CI$trial_num, rev(CI$trial_num)),
             y = c( lo, rev(hi)),
             border=NA,
             col = colors$tr[group_no])
    
  }
  
  
  
  for (group_no in c(1:length(groups))) {
    
    group <- groups[group_no]
    
    
    df <- read.csv(sprintf('data/%s/%s_training_reachdevs.csv', group, group),
                   stringsAsFactors = F)
    
    
    avg <- aggregate(reachdeviation_deg ~ trial_num,
                     data = df,
                     FUN = mean,
                     na.rm = T)
    
    lines(avg,
          col = colors$op[group_no])
    
  }
  
  legend( x = 70,
          y = 12,
          legend = groups,
          col = colors$op[c(1:length(groups))],
          bty='n',
          lty=1,
          # title='groups:',
          bg='#FFFFFF')
  
  axis(side = 1,
       at = c(1,30,60,90))
  axis(side = 2,
       at = c(0,10,20,30))
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}


plotFittedExponentials <- function(target='inline') {
  
  setupFigureFile(target=target,
                  width = 8,
                  height=6,
                  dpi=300,
                  sprintf('doc/fig3_exponentials.%s', target))
  
  groups <- c('control', 'cursorjump', 'handview')
  
  colors <- getColors()
  
  plot(NA, NA,
       main='learning', xlab='trial', ylab='reach deviation [°]',
       xlim=c(0,91), ylim=c(-10,40),
       ax=F, bty='n')
  
  lines(x=c(0,91),
        y=c(30,30),
        col='#999999',
        lw=2,
        lty=2)
  lines(x=c(0,91),
        y=c(0,0),
        col='#999999',
        lw=2,
        lty=2)
  
  
  timepoints <- seq(0,89,0.1)
  
  groups <- c('control', 'cursorjump', 'handview')
  
  # loop through groups
  for (group_no in c(1:length(groups))) {
    
    group <- groups[group_no]
    
    # load the bootstrapped exponential fits:
    expfits <- read.csv(sprintf('data/%s/%s_expfits.csv', group, group),
                        stringsAsFactors = FALSE)
    
    # calculate the 95% CI for lambda and N0
    lambdaCI <- quantile(expfits$lambda, c(0.025, 0.975))
    N0CI     <- quantile(expfits$N0, c(0.025, 0.975))
    
    
    N0 <- unname(N0CI[1])
    lambda <- unname(lambdaCI[1])
    par <- c('N0'=N0,
             'lambda'=lambda)
    fittedcurve <- exponentialModel(par=par,
                                    timepoints=timepoints)
    
    X <- timepoints + 1
    Y <- fittedcurve$output
    
    N0 <- unname(N0CI[2])
    lambda <- unname(lambdaCI[2])
    par <- c('N0'=N0,
             'lambda'=lambda)
    fittedcurve <- exponentialModel(par=par,
                                    timepoints=timepoints)
    
    X <- c(X , rev(timepoints)+1)
    Y <- c(Y , rev(fittedcurve$output))
    
    
    polygon( x = X,
             y = Y,
             border=NA,
             col = colors$tr[group_no])
    
  }
  
  
  # loop through groups
  for (group_no in c(1:length(groups))) {
    
    group <- groups[group_no]
    
    # load the bootstrapped exponential fits:
    expfits <- read.csv(sprintf('data/%s/%s_expfits.csv', group, group),
                        stringsAsFactors = FALSE)
    
    # calculate the 95% CI for lambda and N0
    lambdaCI <- quantile(expfits$lambda, c(0.025, 0.50, 0.975))
    N0CI     <- quantile(expfits$N0, c(0.025, 0.50, 0.975))
    
    N0 <- unname(N0CI[2])
    lambda <- unname(lambdaCI[2])
    
    par <- c('N0'=N0,
             'lambda'=lambda)
    
    fittedcurve <- exponentialModel(par=par,
                                    timepoints=timepoints)
    
    
    lines(x=timepoints+1,
          y=fittedcurve$output,
          col = colors$op[group_no])
    
  }
  
  legend( x = 70,
          y = 12,
          legend = groups,
          col = colors$op[c(1:length(groups))],
          bty='n',
          lty=1,
          # title='groups:',
          bg='#FFFFFF')
  
  axis(side = 1,
       at = c(1,30,60,90))
  axis(side = 2,
       at = c(0,10,20,30))
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}


# nocursor plots -----

plotNoCursors <- function(target='inline') {
  
  setupFigureFile(target=target,
                  width = 8,
                  height=6,
                  dpi=300,
                  sprintf('doc/fig4_nocursors.%s', target))
  
  groups <- c('control', 'cursorjump', 'handview')
  
  colors <- getColors()
  
  plot(NA, NA,
       main='no cursors', xlab='strategy', ylab='reach deviation [°]',
       xlim=c(0,3), ylim=c(-10,40),
       ax=F, bty='n')
  
  lines(x=c(0,3),
        y=c(30,30),
        col='#999999',
        lw=2,
        lty=2)
  lines(x=c(0,3),
        y=c(0,0),
        col='#999999',
        lw=2,
        lty=2)
  
  for (group_no in c(1:length(groups))) {
    
    group <- groups[group_no]
    
    
    df <- read.csv(sprintf('data/%s/%s_nocursors_reachdevs.csv', group, group),
                   stringsAsFactors = F)
    
    
    CI <- aggregate(reachdeviation_deg ~ strategy,
                    data = df,
                    FUN = Reach::getConfidenceInterval)
    
    lo <- CI$reachdeviation_deg[,1]
    hi <- CI$reachdeviation_deg[,2]
    
    polygon( x = c( CI$strategy, rev(CI$strategy))+1,
             y = c( lo, rev(hi)),
             border=NA,
             col = colors$tr[group_no]) 
    
    avg <- aggregate(reachdeviation_deg ~ strategy,
                     data = df,
                     FUN = mean,
                     na.rm = T)
    # print(avg)
    lines(x=c(1,2), avg$reachdeviation_deg, 
          col = colors$op[group_no])
    
    
    participant_avg <- aggregate(reachdeviation_deg ~ strategy + participant,
                                 data = df,
                                 FUN = mean,
                                 na.rm = T)
    
    without <- participant_avg$reachdeviation_deg[which(participant_avg$strategy == 0)]
    with <- participant_avg$reachdeviation_deg[which(participant_avg$strategy == 1)]
    
    points(x = rep(group_no/4, length(without)),
           y = without,
           pch=16,
           cex=1.5,
           col = colors$tr[group_no])
    
    points(x = 2+rep(group_no/4, length(with)),
           y = with,
           pch=16,
           cex=1.5,
           col = colors$tr[group_no])
    
  }
  
  axis(side = 1,
       at = c(1,2),
       labels=c('without','with'))
  
  axis(side = 2,
       at = c(0,15,30))
  
  
  legend( x = 0,
          y = 45,
          legend = groups,
          col = colors$op[c(1:length(groups))],
          bty='n',
          lty=1,
          # title='groups:',
          bg='#FFFFFF',
          xpd=TRUE)
  
}
