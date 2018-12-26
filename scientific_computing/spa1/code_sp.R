##code for doing practical 1


####start by defining functions we will need

get_grade = function(x){
  #Given a fraction of correctly answered questions, returns a letter grade
  x = floor(100*x)
  if (x >= 70){return('A')}
  else if ( x >= 60 ){return('B')}
  else if ( x >= 50 ){return('C')}
  else if ( x >= 40 ){return('D')}
  else {return('F')}
}

cheat_calc = function( a1, a2, c = crib, pc = -log10(pcorr^2), pf = -log10( ((1-pcorr)/4)^2 )){
  #given penalty scores for identical correct (pc) and incorrect (pf) answers,
  #calculates a similarity score between two students with answer vectors a1 and a2
  corrs = (a1 == a2 & a1 == c)
  falses = (a1 == a2 & a1 != c & a1 != 'NA')
  score = pc * sum(corrs) + pf * sum(falses)
return(score)
}


getASCII = function(x){
  #returns TRUE if word x contains only ASCII characters, FALSE otherwise
  x = utf8ToInt(x)
  return(all(x < 128)) #ASCII characters are the first 127 utf8 characters
}

getScore = function(x, scores = sorted){
  #calculates the scrabble score of word x given an alphabetical list of character scores (scores)
  score = 0
  for (i in 1:nchar(x)){ score = score + scores[substr(x, i, i)] }
  return(score)
}

revcomp = function(x, revcomplist = comp){
  #Given a word x, generates the reverse complement of x. revcomplist is a list of reverse complements
  rc = ''
  for (i in nchar(x):1){rc = paste0(rc, revcomplist[substr(x, i, i)])}
  return(rc)
}

can_make = function(x, ls = c('A', 'F', 'L', 'U', 'Y', 'P', 'L', 'N', 'I') ){
  #Given a word x, returns TRUE if it can be made using characters in ls, FALSE otherwise.
  for (i in 1:nchar(x)){
    N = substr(x, i, i)
    if (N %in% ls){ ls = ls[-match(N, ls)]} #see if character in list of allowed characters
    else {return(FALSE)} #return false if no match
  }
  return(TRUE)
}

catch_cheat = function(scores, pairs, cutoff = 0.01){
  #given a list of scores and pairs, reports cheaters with a probability threshold of cutoff
  #Use: catch_cheat(sim_scores, pairs) after running code_sp.R
  p = 0
  while( p < cutoff ){
    top = scores[ order(scores, decreasing = 1)[1] ]
    ind = 
      pair = pairs[ scores == top ]
    pairs = pairs[ scores != top ]; scores = scores[ scores != top ]
    m = mean(scores)
    sd = sd(scores)
    p = pnorm( (m-top)/sd )*(length(scores)+1)
    if (p < cutoff){
      cat('\nstudents', pair, 'appear to be cheating')
      cat('\nsimilarity is', top, 'corresponding to p =', p, '\n')
    }
    else{cat('\nno more students appear to be cheating, most suspicious is', pair, top, p)}
  }
}

###now do the actual exercises

#load file, convert to uppercase and take only unique elements
f = unique(toupper(as.vector(unlist(read.delim('files/usr-share-dict-words', header=FALSE)))))

cat('there are', length(f), 'unique words\n')

aps = grep("'", f) #find words with apostrophes
cat('there are', length(aps), 'words with apostrophes\n')

f = f[-aps] #remove words with apostrophes

ascii = sapply(f, getASCII) #get vector where each element is 1 if all ASCII, 0 otherwise
cat('there are', length(f[!ascii]), 'words with non-ASCIIs')

f = f[ascii] #remove non-ASCII words
 
og = f[grep('*OG$', f)]
ue = f[grep('*OGUE$', f)]

ogues = vector()
for (word in ue){
  #for each XOGUE word, seach our list of YOG words for X=Y
  cat('\n', word, substr(word, 1, nchar(word)-2), '\n')
  w = og[ og == substr(word, 1, nchar(word)-2) ]
  cat(w, '\n')
  ogues = c(ogues, w)
}

cat('\nThere are', length(ogues), 'ogue words. These are\n', ogues, '\n')

scores = read.table('files/scrabble.txt', header=FALSE)
alph = unlist(scores[,1]) #extract list of letters in table
sorted = numeric(25); sorted[as.integer(alph)] = unlist(scores[,4]) #factors are alphabetical
names(sorted) = sort(alph) #now alphabetical list of letter (name) and score (value)

scrabbles = sapply(f, getScore) #compute vector of scrabble scores

png(file='scrabble_hist.png', w=600, h=600)
hist(scrabbles, xlab = 'Scrabble Score', main = 'Score Distribution') #plot hist
dev.off()

m = max(scrabbles) #max score
mword = names(scrabbles[which(scrabbles == m)]) #highest scoring word
cat('\nmax word is', mword, 'with', m,'points')

comp = rev(as.vector(sort(alph))); names(comp) = rev(comp) #vector of reverse complement letters
revcomps = vector()
#find all words with a reverse complement in the database
for (w in f){ if (revcomp(w) %in% f){revcomps = c(revcomps, w)}}

temp = revcomps
revcomps = sapply(temp, revcomp)
#generate list where names are words, values are reverse complements of those words
names(revcomps) = temp
cat('\n\n', revcomps, '\n')


####Now find the possible word combinations.

#First construct new database of allowed words, then check if they can be made from our characters

witha = f[ grep('A', f) ]; witha = witha[ nchar(witha) > 3 ]
ls = c('F', 'A', 'L', 'U', 'Y', 'P', 'L', 'N', 'I') #allowed characters
#test = c('APLN', 'ABBA', 'AFLUY') #test case

realwords = witha[ as.vector(unlist(sapply(witha, can_make))) ] #get possible words
cat('\nwe can construct', length(realwords), 'unique words from our letters. These are:\n')
print(realwords)


####Now do the student marking and cheating task

#read in data
crib = as.vector(unlist(read.delim('files/grading/crib.dat', header=FALSE)))
grade_scheme = read.table('files/grading/grade.txt', header=TRUE)

#initialize dataframe of results
res = data.frame(student = 1:12, score = numeric(12), grade = character(12), rank = numeric(12),
                  stringsAsFactors = FALSE)
answers = vector('list', 12) #list of vectors of answers
scores = numeric(12) #vector of student scores
for (i in 1:12){ #treat each student individually
  ans = read.table(paste0('files/grading/student',i,'.dat'), header=TRUE, stringsAsFactors = FALSE)
  #make a list of 100 'NA', change the answered questions to the relevant answer
  answers[[i]] = rep('NA', 100); answers[[i]][ans[,'qn']]=ans[,'response']
  score = 0
  for (j in 1:dim(ans)[1]){ if (ans[j,2] == crib[ans[j,1]]){score = score+1} } #count scores
  cat('\n', i, 'score', score)
  res[i, 'score' ] = score; scores[i] = score
  cat('\n', 'grade', get_grade(score/30), '\n')
  res[i, 'grade' ] = get_grade(score/30)
}

scores = sort(scores, decreasing = TRUE) #sort by score to find ranks
for (i in 1:12){res[i, 'rank'] = which(scores == res[i, 'score'])[1]} #set ranks

###Test for cheating

pcorr = mean(scores)/30 #total probability of a question being answered correctly
pc = -log10(pcorr^2) #probabiliy that two students both provide a correct answer
pf = -log10( ((1-pcorr)/4)^2 ) #probability that two students both provide the same incorrect answer

cheats =  matrix(NA, 12, 12)  #initialize matrix of similarity scores
plotCheat = matrix(0, 12, 12) #need a slightly different data structure for plotting using 'image'
sim_scores = numeric(0); pairs = numeric(0) #construct these for automatic cheat detection

for (i in 1:11){
  for (j in (i+1):12){
    cheat = cheat_calc( answers[[i]], answers[[j]]) #calculate pairwise similarity scores
    cheats[i,j] = cheats[j,i] = cheat
    plotCheat[13-i, j] = plotCheat [13-j, i] = cheat
    sim_scores = c(sim_scores, cheat); pairs = c(pairs, paste0('(',i,',',j,')') )
  }
}


png(file='cheat_heat.png', w=600, h=600)
image(t(plotCheat), col=heat.colors(20), axes = 0) #plot heatmap
axis(3, at = seq(0, 1, length = 12), labels=1:12,srt=45,tick=FALSE)
axis(2, at = seq(0, 1, length = 12), labels=12:1,srt=45,tick=FALSE)
dev.off()

u = sim_scores[ order(sim_scores, decreasing = 1)[1:5] ] #find the top similarity scores
cat( u )

inds = vector('list', length(u))
for (i in 1:length(u) ){ inds[[i]] = which( cheats == u[i], arr.ind = 1 )[1,] }
print(inds) #find the most similar pairs

png(file='cheat_hist.png', w=600, h=600) #plot histogram.
hist(sim_scores, 30, main = 'Distribution of Similarity Scores',
     xlab = 'Similarity', ylab = 'Frequency')
dev.off()

sds = numeric(5)
for (i in 1:5){
  #get the number of standard deviations of the top similarity scores from the mean of the similarity
  #distribution excluding that particular score
  m = mean(sim_scores[sim_scores != u[i]], na.rm = 1)
  sd = sd(sim_scores[sim_scores != u[i]], na.rm = 1)
  sds[i] = (u[i]-m)/sd
  cat('\n', m, sd, sds[i], '\n')
}

new_scores = sim_scores[sim_scores != u[[1]]] #remove highest scoring score
for (i in 2:5){
  #get the number of standard deviations of the top similarity scores from the mean of the similarity
  #distribution excluding that particular score after removing 2,6
  m = mean(new_scores[new_scores != u[i]], na.rm = 1)
  sd = sd(new_scores[new_scores != u[i]], na.rm = 1)
  sds[i] = (u[i]-m)/sd
  cat('\n', m, sd, sds[i], '\n')
}


#####Or we can do the above 'automatically', listing only the cheating results

catch_cheat(sim_scores, pairs)
   
   