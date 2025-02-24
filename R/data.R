# download data from OSF -----

getData <- function() {
  
  Reach::downloadOSFdata(repository  = 'm5dt4',
                         filelist    = list('data/' = c(
                           'control_aligned.zip',
                           'control_rotated.zip',
                           'handview_aligned.zip',
                           'handview_rotated.zip',
                           'cursorjump_aligned.zip',
                           'cursorjump_rotated.zip',
                           'demographics.csv',
                           'schedule.csv'
                         )
                         ),
                         folder      = 'data/',
                         overwrite   = TRUE,
                         unzip       = TRUE,
                         removezips  = TRUE)
  
  Reach::downloadOSFdata(repository  = 'm5dt4',
                         filelist    = list(
                           'documents/' = c('explicit_task_order.png') 
                         ),
                         folder      = 'doc/',
                         overwrite   = TRUE,
                         unzip       = TRUE,
                         removezips  = TRUE)
}


# participants -----

# get a list of participants that were in a particular group / condition:

groupParticipants <- function(group) {
  
  demo <- read.csv('data/demographics.csv', stringsAsFactors = F)
  
  participants <- demo$participant[which(demo$group == group)]
  
  return(participants)
  
}

# schedule -----

getSchedule <- function(order=NULL, addsubtasktrialnums=FALSE) {
  
  schedule <- read.csv('data/schedule.csv', stringsAsFactors = F)
  
  # figure out strategy include/exclude order, and replace strategy column with 
  # 'include' and 'exclude' values
  if (!is.null(order)) {
    if (is.character(order) & length(order) == 1) {
      if (order %in% c('EI', 'IE')) {
        # valid order
        
        if (order == 'IE') {
          schedule$strategy[which(schedule$strategy == 0)] <- 'include'
          schedule$strategy[which(schedule$strategy == 1)] <- 'exclude'
        }
        if (order == 'EI') {
          schedule$strategy[which(schedule$strategy == 0)] <- 'exclude'
          schedule$strategy[which(schedule$strategy == 1)] <- 'include'
        }
        
        
      }
    }
  }
  # when order is not provided (correctly) the strategy values will not be set
  
  # if the callers wants subtask trial numbers as well, they can be added here:
  if (addsubtasktrialnums) {
    
    schedule$subtask_trial <- 0
    
    for (subtask_no in unique(schedule$subtask)) {
      
      idx <- which(schedule$subtask == subtask_no)
      
      schedule$subtask_trial[idx] <- c(1:length(idx))
      
    }
    
  }
  
  # return the current schedule
  
  return(schedule)
  
}

