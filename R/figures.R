# general plotting functions----
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

#actual plotting ----
plotTraining <- function(target='inline'){
  
  setupFigureFile(target=target,
                  width = 8,
                  height = 6,
                  dpi = 300,
                  sprintf('doc/fig2_training.%s', target))
  
  groups <- c('control', 'cursorjump', 'handview')
  colors <- getColors()
  
  plot(NA, NA,
       main = 'learning', xlab = 'trial', ylab = 'reach deviation (deg)',
       xlim = c(0, 91), ylim = c(-10, 40),
       ax = F, bty = 'n')
  
  lines(x=c(0,91), y=c(30,30), col='#999999', lty=2, lwd=2)
  
  for(group_no in c(1:length(groups))){
    
    group <- groups[group_no]
    
    df <- read.csv(filename <- sprintf('data/%s/%s_training_reachdevs.csv', group, group), stringsAsFactors = F)
    
    CI <- aggregate(reachdeviation_deg ~ trial_num,
                     data = df, 
                     FUN = Reach::getConfidenceInterval)
    
    lo <- CI$reachdeviation_deg[,1]
    hi <- CI$reachdeviation_deg[,2]
    
    polygon(x = c(CI$trial_num, rev(CI$trial_num)),
            y = c(lo, rev(hi)),
            border = NA,
            col = colors$tr[group_no])
    
  }
  
  for(group_no in c(1:length(groups))){
    
    group <- groups[group_no]
    
    df <- read.csv(filename <- sprintf('data/%s/%s_training_reachdevs.csv', group, group), stringsAsFactors = F)
    
    avg <- aggregate(reachdeviation_deg ~ trial_num,
                     data = df, 
                     FUN = mean,
                     na.rm = T)
    
    lines(avg,
          col = colors$op[group_no])
    
  }
  
  legend(x = 60, y = 15,
         legend = groups,
         col = colors$op[c(1:length(groups))],
         bty = 'n',
         lty = 1, lwd = 2)#,
         #title = 'groups:')
  
  axis(side = 1,
       at = c(1, 30, 60, 90))
  axis(side = 2,
       at = c(0, 10, 20, 30))
  
  if(target %in% c('pdf', 'svg', 'png', ' tiff')){
    dev.off()
  }
  
}