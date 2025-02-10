source('R/data.R')

getParticipantTraining <- function(group, participant){
  
  AL_file <- sprintf('data/%s/%s/%s_aligned_training.csv', group, participant, participant)
  RO_file <- sprintf('data/%s/%s/%s_rotated_training.csv', group, participant, participant)
  
  AL_df <- read.csv(file = AL_file,
                    stringsAsFactors = F)
  RO_df <- read.csv(file = RO_file,
                    stringsAsFactors = F)
  
  return(list('aligned' = AL_df,
         'rotated' = RO_df))
  
}

getGroupTraining <- function(group){
  
  participants <- groupParticipants(group = group)
  
  for (participant in participants[c(1:3)]){
    
    participant_df <- getParticipantTraining(group = group, participant = participant)
    
    baseline <- getBaseline(participant_df[['aligned']])
    
  }
  
}

getBaseline <- function(df){
  
  
  df <- df[which(df$trial_num %in% c(31:45)),]
  str(df)
  
  trialnos <- unique(df$trial_num)
  
  for(trial in trialnos){
    
    tdf <- df[which(df$trial_num == trial),]
    reachdev <- getReachDeviation(tdf)
    print(reachdev)
    # 
    # target <- tdf$targetangle_deg[1]
    # 
    # X <- tdf$handx_cm
    # Y <- tdf$handy_cm
    
    
    
  }
  
  
  
}

getReachDeviation <- function(df){
  
  target <- df$targetangle_deg[1]
  
  X <- tdf$handx_cm
  Y <- tdf$handy_cm
  
  # at 1/4 target distance
  #target distance was 10 cm
  #2.5 cm
  distances <- sqrt(X^2 + Y^2)
  idx <- which(distances > 2.5)[1]
  
  x <- X[idx]
  y <- Y[idx]
  
  th <- (-1*target/180) * pi
  R <- t(matrix(data=c(cos(th),sin(th),-sin(th),cos(th)),nrow=2,ncol=2))
  
  
  # rotate the coordinates, add the origin back in
  norm_sample <- matrix(data=c(x,y), ncol=2) %*% R
  
  #print(norm_sample)
  
  reachdev <- (atan2(norm_sample[2], norm_sample[1]) / pi) *180
  #print(reachdev)
  return(c('trial_num' = df$trial_num[1], 'targetangle_deg' = target, 'reachdeviation_deg' = reachdev))
}

















