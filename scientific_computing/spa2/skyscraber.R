###Code for sp2; generate knitr afterwards.

skyplot = function(edges, soln = '', ...){
  
  N = length(edges) / 4
  edges[edges == 0] = ''; soln[soln==0]='' #zero values indicate blanks
  s = 0.8/N #step length
  cex = 7/N
  plot.new()
  title(...)
  segments( c(0.1,0.1,0.1,0.9), c(0.1,0.9,0.1,0.1), c(0.9,0.9,0.1,0.9), c(0.1,0.9,0.9,0.9) ) #outline
  segments(0.1, 0.1+(1:N)*s, 0.9, 0.1+(1:N)*s) #horizontal lines
  segments(0.1+(1:N)*s, 0.1, 0.1+(1:N)*s, 0.9) #vertical lines
  
  pos = c(3, 4, 1, 2)
  #add top + bottom text
  text( s*(0:(N-1)+0.5)+0.1, y=rep(c(.91, 0.09),each=N),
        edges[c(1:N,(3*N):(2*N+1))], pos = rep(c(3,1),each=N), cex = cex )
  #add left+right text
  text( rep(c(0.91,0.09),each=N), y=s*(0:(N-1)+0.5)+0.1,
        edges[c((2*N):(N+1),(3*N+1):(4*N))], pos = rep(c(4,2),each=N), cex = cex )
  text( rep(s*(0:(N-1)+0.5)+0.1, N), y=rep(s*((N-1):0+0.5)+0.1, each=N),
        soln, cex = cex)
}

permutations <- function(n){
  #returns a matrix or all possible permutations of the numbers 1:n
  if(n==1){
  return(matrix(1))
    } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

visible = function(vec, reverse = 0){
  #given a vector of skyscraper heights, returns the number visible ones
  #if reverse, reverses the order of the skyscrapers
  if (reverse) vec = rev(vec) #look from other direction
  n = 0; max = 0
  for (sky in vec){ if (sky > max){ n = n+1; max = sky } }
  return(n)
}

all_visible = function(mat){
  #given a matrix, generates a vector of the number of visible skyscapers
  vis = 
    c(apply(mat, 2, visible), #consider top edges
    apply(mat, 1, visible, reverse = 1), #right edges
    rev(apply(mat, 2, visible, reverse = 1)), #bottom edges, reversing output order
    rev(apply(mat, 1, visible))) #left edges, reversing output order
  return(vis)
}

solve_sky = function(edges,main = 'puzzle'){
  #given a list of edge values of visible skyscrapers for a 4x4 problem,
  #finds, prints and plots a solution to the problem
  sol = matrix(,4,4)
  sols = list()
  perms1 = permutations(4)
  new_main = main
  #print(perms1)
  
  for (i in 1:factorial(4)){ #consider all permutations of 1:4
    sol[1,] = perms1[i,]
    #at most 6 solutions
    #we can now remove the permurations that are not allowed together with
    #the first permutation (gives Nsol = 24*9*4*1 = 864 rather than 24^4 = 331776)
    perms2 = perms1[apply(perms1, 1, function(x){ !any(x == sol[1,]) }), ]
    
    for (j in 1:dim(perms2)[1]){
      sol[2,] = perms2[j,] #add next row
      perms3 = perms2[apply(perms2, 1, function(x){ !any(x == sol[2,]) }), ]
      
      for (k in 1:dim(perms3)[1]){
        sol[3,] = perms3[k,] #add next row
        #now we can infer final row
        sol[4,] = perms3[apply(perms3, 1, function(x){ !any(x == sol[3,]) }), ]
        #print(sol)
        vis = all_visible(sol) #get list of visible skyscrapers for sol
        #print(vis)
        vis[ edges == 0 ]=0 #don't look at edges with no information
        
        if ( all(edges == vis) ){ #check if solution consistent with problem
          cat('\nfound solution\n')
          print(sol)
          L = length(sols)
          if (L>0) new_main = paste0(main,as.character(L+1))
          skyplot(edges,c(t(sol)),main=new_main)
          sols[[L+1]] = c(t(sol))
        }
      }
    }
  }
  N = length(sols)
  if (N > 0){
    cat('\nfound', N, 'solutions\n\n')
    return(sols)
  }else{cat('\npuzzle not solvable, sorry\n')} #if our proposed problem is unsolvable
}

#for (N in 4:9) skyplot(1:(4*N), 1:(N^2))

par(mfrow=c(1,2),mar=c(2,2,1.5,1))
skyplot( 1:20, 1:25, main = 'A')
skyplot( c(2,0,0,0,2,  0,0,0,2,0,  0,0,3,4,0,  0,1,0,2,0),
         c(0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,5,  5,0,0,0,4, 0,0,0,0,0),
         main='B')

par(mfrow=c(1,2),mar=c(2,2,1.5,1))

e1 = scan(comment.char="#", quiet=T,
          'https://raw.githubusercontent.com/sje30/rpc2018/master/a2/e1.dat')
e2 = c(3,0,1,0, 0,0,1,0, 0,0,0,0, 0,0,0,2)

a = solve_sky(e1, main = 'A')
b = solve_sky(e2, main = 'B')

testm = c(0,2,1,0, 0,0,1,3, 0,2,0,0, 1,0,2,2)

c = solve_sky(testm, main = 'test')
