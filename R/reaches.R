getParticipantTraining <- function(group, participant){
  
  AL_file <- sprintf('data/%s/%s/%s_aligned_training.csv', group, participant, participant)
  RO_file <- sprintf('data/%s/%s/%s_rotated_training.csv', group, participant, participant)
  
  AL_df <- read.csv(file = AL_file,
                    stringsAsFactors = F)
  RO_df <- read.csv(file = RO_file,
                    stringsAsFactors = F)
  
}

getGroupTraining <- function(group){
  
  participants <- groupParticipants(group = group)
  
  for (participant in participants){
    
    participant_df <- getParticipantTraining(group = group, participant = participant)
    
  }
  
}