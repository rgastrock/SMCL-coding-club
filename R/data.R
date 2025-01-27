

# download data from OSF:
getData <- function(){
  
  Reach::downloadOSFdata(repository = 'm5dt4',
                         filelist = list('data/' =c('control_aligned.zip', 'control_rotated.zip', 'handview_aligned.zip', 'handview_rotated.zip', 'cursorjump_aligned.zip', 'cursorjump_rotated.zip')),
                         folder = 'data/',
                         overwrite = TRUE,
                         unzip = TRUE,
                         removezips = TRUE)
  
  
}