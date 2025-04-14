

# training reach deviations ----

getAllTraining <- function() {
  
  groups <- c('control', 'cursorjump', 'handview')
  
  for (group in groups) {
    
    df <- getGroupTraining(group)
    
    filename <- sprintf('data/%s/%s_training_reachdevs.csv', group, group)
    
    write.csv( df,
               filename,
               row.names=FALSE)
    
  }
  
}



getParticipantTraining <- function(group, participant) {
  
  AL_file <- sprintf('data/%s/%s/%s_aligned_training.csv', group, participant, participant)
  RO_file <- sprintf('data/%s/%s/%s_rotated_training.csv', group, participant, participant)
  
  AL_df <- read.csv( file = AL_file,
                     stringsAsFactors = F)
  RO_df <- read.csv( file = RO_file,
                     stringsAsFactors = F)
  
  return( list('aligned'=AL_df,
               'rotated'=RO_df) )
  
}


getGroupTraining <- function(group) {
  
  participants <- groupParticipants(group = group)
  
  # print(participants)
  
  learningCurves <- NA
  
  for (participant in participants) {
    
    participant_df <- getParticipantTraining( group       = group,
                                              participant = participant )
    
    # print(participant_df)
    
    baseline <- getBaseline( df = participant_df[['aligned']] )
    
    baseline <- removeOutliers(baseline, rotation = 0)
    
    baseline <- aggregate(reachdeviation_deg ~ targetangle_deg, data=baseline, FUN=median, na.rm=TRUE)
    
    rotated <- getRotatedLearning( df = participant_df[['rotated']] )
    
    rotated <- removeOutliers(rotated, rotation = -30)
    
    rotated <- baselineCorrection(baseline=baseline, rotated=rotated)
    
    rotated$participant <- participant
    
    if (is.data.frame(learningCurves)) {
      learningCurves <- rbind(learningCurves, rotated)
    } else {
      learningCurves <- rotated
    }
    
    # str(rotated)
    
    # print(participant)
    # print(length(which(is.na(rotated$reachdeviation_deg))))
    # 
    # if (participant %in% c("3c0021","c831b7")){
    #   plot(rotated$reachdeviation_deg)
    # }
    
    
    # rotated <- getRotatedLearning( df = participant_df[['rotated']] )
    
    # removeOutliers(rotated, rotation=30)
    
  }
  
  # plot(x=learningCurves$trial_num,
  #      y=learningCurves$reachdeviation_deg)
  
  return(learningCurves)
  
}

getBaseline <- function(df, task='training') {
  
  schedule <- read.csv('data/schedule.csv', stringsAsFactors = F)
  subtasks <- unique(schedule[which(schedule$session == 'aligned' & schedule$task == 'training'),]$subtask)
  
  if (task == 'training') {
    trialnums <- c(31:45)
    for (subtask in subtasks[c(2:length(subtasks))]) {
      sttn <- schedule$trial_num[which(schedule$subtask == subtask)]
      trialnums <- c(trialnums, sttn[7:9])
    }
  }
  if (task == 'nocursor') {
    trialnums <- unique(df$trial_num)
    # trialnums <- trialnums[which(trialnums < 108)]
  }
  
  # str(df)
  df <- df[which(df$trial_num %in% trialnums),]
  # str(df)
  
  
  
  trialnos <- unique(df$trial_num)
  
  outdf <- NA
  
  for (trial in trialnos) {
    
    tdf <- df[which(df$trial_num == trial),]
    
    reachdev <- getReachDeviation(tdf)
    
    reachdev <- data.frame(t(data.frame(reachdev)))
    
    if (is.data.frame(outdf)) {
      outdf <- rbind(outdf, reachdev)
    } else {
      outdf <- reachdev
    }
    
  }
  
  return(outdf)
  
}


getRotatedLearning <- function(df, task='training') {
  
  if (task == 'training') {
    trialnums <- c(1:90)
  }
  
  if (task == 'nocursor') {
    trialnums <- unique(df$trial_num)
    # trialnums <- trialnums[which(trialnums < 108)]
  }
  
  df <- df[which(df$trial_num %in% trialnums),]
  
  trialnos <- unique(df$trial_num)
  
  outdf <- NA
  
  for (trial in trialnos) {
    
    tdf <- df[which(df$trial_num == trial),]
    
    reachdev <- getReachDeviation(tdf)
    
    reachdev <- data.frame(t(data.frame(reachdev)))
    
    if (is.data.frame(outdf)) {
      outdf <- rbind(outdf, reachdev)
    } else {
      outdf <- reachdev
    }
    
  }
  
  return(outdf)
  
}


getReachDeviation <- function(df) {
  
  
  target <- df$targetangle_deg[1]
  
  X <- df$handx_cm
  Y <- df$handy_cm
  
  # at 1/4 the target distance
  # target distance was 10 cm
  # 2.5 cm
  
  distances <- sqrt(X^2 + Y^2)
  idx <- which(distances > 2.5)[1]
  
  # print(idx)
  x <- X[idx]
  y <- Y[idx]
  
  th <- (-1*target/180) * pi
  R <- t(matrix(data=c(cos(th),sin(th),-sin(th),cos(th)),nrow=2,ncol=2))
  
  # rotate the coordinates, add the origin back in
  norm_sample <- matrix(data=c(x,y),ncol=2) %*% R
  
  # print(norm_sample)
  
  reachdev <- (atan2(norm_sample[2], norm_sample[1]) / pi) * 180
  
  # print(reachdev)
  
  return(c('trial_num' = df$trial_num[1],
           'targetangle_deg'=target,
           'reachdeviation_deg'=reachdev))
  
}

removeOutliers <- function(df, rotation=0) {
  
  windowwidth <- 50
  
  if (rotation == 0) {
    
    df$reachdeviation_deg[which(abs(df$reachdeviation_deg) > windowwidth)] <- NA
    # df$reachdeviation_deg[which(df$reachdeviation_deg >  windowwidth)] <- NA
    # df$reachdeviation_deg[which(df$reachdeviation_deg < -windowwidth)] <- NA
    
  } else {
    
    hi <-  windowwidth
    lo <- -windowwidth
    
    if (rotation > 0) {
      lo <- lo - rotation
    } else {
      hi <- hi - rotation
    }
    
    df$reachdeviation_deg[which(df$reachdeviation_deg > hi)] <- NA
    df$reachdeviation_deg[which(df$reachdeviation_deg < lo)] <- NA
    
  }
  
  return(df)
  
}


baselineCorrection <- function(baseline=baseline, rotated=rotated) {
  
  for (target in baseline$targetangle_deg) {
    
    bias <- baseline$reachdeviation_deg[which(baseline$targetangle_deg == target)]
    idx <- which(rotated$targetangle_deg == target)
    rotated$reachdeviation_deg[idx] <- rotated$reachdeviation_deg[idx] - bias
    
  }
  
  return(rotated)
  
}


# no cursor reach deviations ----


getAllNoCursors <- function() {
  
  groups <- c('control', 'cursorjump', 'handview')
  
  for (group in groups) {
    
    df <- getGroupNoCursors(group)
    
    filename <- sprintf('data/%s/%s_nocursors_reachdevs.csv', group, group)
    
    write.csv( df,
               filename,
               row.names=FALSE)
    
  }
  
}


getParticipantNoCursors <- function(group, participant) {
  
  AL_file <- sprintf('data/%s/%s/%s_aligned_nocursor.csv', group, participant, participant)
  RO_file <- sprintf('data/%s/%s/%s_rotated_nocursor.csv', group, participant, participant)
  
  AL_df <- read.csv( file = AL_file,
                     stringsAsFactors = F)
  RO_df <- read.csv( file = RO_file,
                     stringsAsFactors = F)
  
  return( list('aligned'=AL_df,
               'rotated'=RO_df) )
  
}


getGroupNoCursors <- function(group) {
  
  participants <- groupParticipants(group = group)
  
  NoCursorReaches <- NA
  
  for (participant in participants) {
    
    participant_df <- getParticipantNoCursors( group       = group,
                                               participant = participant )
    
    
    baseline <- getBaseline( df = participant_df[['aligned']],
                             task = 'nocursor')
    
    
    baseline <- removeOutliers(baseline, rotation = 0)
    
    baseline <- aggregate(reachdeviation_deg ~ targetangle_deg, data=baseline, FUN=median, na.rm=TRUE)
    
    
    withStrategy <- participant_df[['rotated']][which(participant_df[['rotated']]$strategy == 1),]
    withoutStrategy <- participant_df[['rotated']][which(participant_df[['rotated']]$strategy == 0),]
    
    withStrategy <- getRotatedLearning( df = withStrategy,
                                        task = 'nocursor')
    withoutStrategy <- getRotatedLearning( df = withoutStrategy,
                                           task = 'nocursor')
    
    withStrategy <- removeOutliers(withStrategy, rotation = -30)
    withoutStrategy <- removeOutliers(withoutStrategy, rotation = -30)
    
    withStrategy <- baselineCorrection(baseline=baseline, rotated=withStrategy)
    withoutStrategy <- baselineCorrection(baseline=baseline, rotated=withoutStrategy)
    
    withStrategy$participant <- participant
    withoutStrategy$participant <- participant
    
    withStrategy$strategy <- 1
    withoutStrategy$strategy <- 0
    
    df <- rbind(withStrategy, withoutStrategy)
    
    if (is.data.frame(NoCursorReaches)) {
      NoCursorReaches <- rbind(NoCursorReaches, df)
    } else {
      NoCursorReaches <- df
    }
    
  }
  
  return(NoCursorReaches)
  
}



getNoCursorData <- function() {
  
  groups <- c('control', 'cursorjump', 'handview')
  
  all_data <- NA
  
  for (group_no in c(1:length(groups))) {
    
    group <- groups[group_no]
    
    df <- read.csv(sprintf('data/%s/%s_nocursors_reachdevs.csv', group, group),
                   stringsAsFactors = F)
    
    participant_avg <- aggregate(reachdeviation_deg ~ strategy + participant,
                                 data = df,
                                 FUN = mean,
                                 na.rm = T)
    
    participant_avg$group <- group
    
    if (is.data.frame(all_data)) {
      all_data <- rbind(all_data, participant_avg)
    } else {
      all_data <- participant_avg
    }
    
  }
  
  return(all_data)
  
}
