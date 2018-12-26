play_A = function(){
  if( runif(1) < (18/37) ){return( c(1, 1) ) #if red, win $1. Always place 1 bet
  }else{return( c(-1,1) )} #lose $1
}

play_B = function(){
  if( runif(1) < (1/37) ){return( c(35, 1) ) #if number, win $35. Always place 1 bet
  }else{return(c(-1,1))} #lose $1
}

play_C = function(show = 0){
  wins = 0; bet = 1; n = 0
  
  while ( wins < 10 & bet <= 100 ){ #continue playing until won 10 or bet > 100
    if(show) cat('\nnet bet', bet, '\n')
    n = n+1 #count number of bets
    if( runif(1) < (18/37) ){ wins = wins+bet; bet = 1 #win amount bet, reset bet
    }else{ wins = wins-bet; bet = 2*bet } #lose amount bet, double bet
  }
  return( c(wins, n) )
}

play_D = function(show = 0){
  nums  = 1:4; wins = 0; n = 0; bet = 5
  
  while (length(nums) > 0 & bet <= 100){ #continue playing until empty list or bet > 100
    if(show) cat('\nnew nums', nums, 'bet', bet, '\n')
    n = n+1 #count number of bets
    
    if( runif(1) < (18/37) ){ wins = wins+bet; nums = nums[-c(1,length(nums))]
    #win amount bet, remove first, last num
    
    }else{ wins = wins-bet; nums = c(nums, bet) }
    #lose amount bet, add sum of first, last num to list
    
    if (length(nums)==0){bet = 0 #we're done playing
    }else if (length(nums) == 1){bet = nums #if only one num, bet that next
    }else{bet = nums[1]+nums[length(nums)]} #bet sum of first and last number}
    
  }
    
  return( c(wins, n) )
}

play_many = function(type, N=100000){
  #Takes as input a type of game (A,B,C,D), plays it N times
  #and returns a list of mean and sd of amount won/lost, games won and bets placed
  f = switch(type, #pick correct game
         'A' = play_A,
         'B' = play_B,
         'C' = play_C,
         'D' = play_D,
         )
  
  amounts = bets = wins = numeric(N)
  for (i in 1:N){
    out = f() #play game and add results to lists
    amounts[i]=out[1]
    bets[i]=out[2]
    wins[i] = as.numeric( out[1] > 0 )
  }
  cat('\nplayed game', type, 'won', 100*mean(wins), '% of games expecting',
      mean(amounts), 'money in', mean(bets), 'bets\n')
  return(c(mean(amounts), mean(wins), mean(bets), sd(amounts), sd(wins), sd(bets)))
}

sim_games = function(){
  #for each of games (A,B,C,D), plays it 100000 times and returns a df with results
  awin = pwin = nbets = numeric(4)
  types = c('A', 'B', 'C', 'D')
  for (type in 1:4){
    out = play_many(types[type]) #play 100000 games and store results
    awin[type] = paste0(as.character(round(out[1],3)),', ', as.character(round(out[4],2)))
    pwin[type] = paste0(as.character(round(100*out[2],1)),'%, ', as.character(round(100*out[5],1)), '%')
    nbets[type] = paste0(as.character(round(out[3],1)),', ', as.character(round(out[6],1)))
  }
  df = data.frame('Winnings_(mean,sd)' = awin,
                  'Prop.wins_(mean,sd)' = pwin,
                  'Plat.time_(mean,sd)' = nbets)
  row.names(df) = types
  return(df)
}

df = sim_games()

repeat_sim = function(){
  #run 100000 game simulations five times and return data frame
  #with min and max of amounts won, proportion won and bets placed.
  wins = props = bets = matrix(,4,5)
  for (i in 1:5){
    df = sim_games()
    wins[,i] = df[,1]; props[,i] = df[,2]; bets[,i] = df[,3] 
  }
  df_summ = data.frame('Win min' = apply(wins, 1, min),
                  'Win max' = apply(wins, 1, max),
                  'Prop min' = apply(props, 1, min),
                  'Prop max' = apply(props, 1, max),
                  'Exp min' = apply(bets, 1, min),
                  'Exp max' = apply(bets, 1, max))
  row.names(df_summ) = row.names(df_summ)
  return(df_summ)
}

#df = repeat_sim()


