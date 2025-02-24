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
  
  learningCurves <- NA
  
  for (participant in participants){
    
    participant_df <- getParticipantTraining(group = group, participant = participant)
    
    baseline <- getBaseline(participant_df[['aligned']])
    
    baseline <- removeOutliers(baseline, rotation = 0)
    baseline <- aggregate(reachdeviation_deg ~ targetangle_deg, data = baseline, FUN = median, na.rm = T)
    
    rotated <- getRotatedLearning (df = participant_df[['rotated']])
    rotated <- removeOutliers(rotated, rotation = -30)
    
    rotated <- baselineCorrection(baseline = baseline, rotated = rotated)
    
    rotated$participant <- participant
    
    if(is.data.frame(learningCurves)){
      learningCurves <- rbind(learningCurves, rotated)
    } else {
      learningCurves <- rotated
    }
    
    # print(baseline$reachdeviation_deg)
  }
  # plot(x=learningCurves$trial_num,
  #      y=learningCurves$reachdeviation_deg)
  return(learningCurves)
}

getBaseline <- function(df){
  
  schedule <- read.csv('data/schedule.csv', stringsAsFactors = F)
  subtasks <- unique(schedule[which(schedule$session == 'aligned' & schedule$task == 'training'),]$subtask)
  trialnums <- c(31:45)
  
  for(subtask in subtasks[c(2:length(subtasks))]){
    sttn <- schedule$trial_num[which(schedule$subtasks == subtask)]
    trialnums <- c(trialnums, sttn[7:9])
  }
  
  df <- df[which(df$trial_num %in% trialnums),]
  # str(df)
  
  trialnos <- unique(df$trial_num)
  
  outdf <- NA
  
  for(trial in trialnos){
    
    tdf <- df[which(df$trial_num == trial),]
    reachdev <- getReachDeviation(tdf)
    reachdev <- data.frame(t(data.frame(reachdev)))
    
    if(is.data.frame(outdf)){
      outdf <- rbind(outdf, reachdev)
    } else {
      outdf <- reachdev
    }
    
    
    
    
  }
  return(outdf)
  
  
}

getRotatedLearning <- function(df){
  
  trialnums <- c(1:90)
  
  df <- df[which(df$trial_num %in% trialnums),]
  # str(df)
  
  trialnos <- unique(df$trial_num)
  
  outdf <- NA
  
  for(trial in trialnos){
    
    tdf <- df[which(df$trial_num == trial),]
    reachdev <- getReachDeviation(tdf)
    reachdev <- data.frame(t(data.frame(reachdev)))
    
    if(is.data.frame(outdf)){
      outdf <- rbind(outdf, reachdev)
    } else {
      outdf <- reachdev
    }
    
    
    
    
  }
  return(outdf)
  
  
}

getReachDeviation <- function(df){
  
  target <- df$targetangle_deg[1]
  
  X <- df$handx_cm
  Y <- df$handy_cm
  
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

removeOutliers <- function(df, rotation = 0){
  
  windowwidth <- 50
  
  if(rotation == 0){
    df$reachdeviation_deg[which(abs(df$reachdeviation_deg) > windowwidth)] <- NA
  } else{
    hi <- windowwidth
    lo <- -windowwidth
    if(rotation > 0){
      lo <- lo - rotation
    } else {
      hi <- hi - rotation
    }
    df$reachdeviation_deg[which(df$reachdeviation_deg > hi)] <- NA
    df$reachdeviation_deg[which(df$reachdeviation_deg < lo)] <- NA
  }
  
  
  return(df)
}

baselineCorrection <- function(baseline = baseline, rotated = rotated){
  
  for(target in baseline$targetangle_deg){
    bias <- baseline$reachdeviation_deg[which(baseline$targetangle_deg == target)]
    idx <- which(rotated$targetangle_deg == target)
    rotated$reachdeviation_deg[idx] <- rotated$reachdeviation_deg[idx] - bias
  }
  
  return(rotated)
  
}













