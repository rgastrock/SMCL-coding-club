

# download data from OSF----
getData <- function(){
  
  Reach::downloadOSFdata(repository = 'm5dt4',
                         filelist = list('data/' =c('control_aligned.zip', 'control_rotated.zip', 'handview_aligned.zip', 'handview_rotated.zip', 'cursorjump_aligned.zip', 'cursorjump_rotated.zip','demographics.csv')),
                         folder = 'data/',
                         overwrite = TRUE,
                         unzip = TRUE,
                         removezips = TRUE)
  
  
}


# participants----

#get a list of participants that were in a particular group/ condition

groupParticipants <- function(group){
  
  demo <- read.csv('data/demographics.csv', stringsAsFactors = F)
  participants <- demo$participant[which(demo$group == group)]
  
  return(participants)
  
}


