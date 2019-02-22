#function implementing dmin2d model
include("plots.jl")

function dmin2d(n, m, s, xlo, xhi, ylo, yhi; hs = false,
                Plot = false, Print = false, thresh = 10000)
  ## n: number of points to simulate
  ## m: mean of Normal distribution
  ## s: s.d. of Normal distribution
  ## xlo, xhi: possible range of X values.
  ## ylo, yhi: possible range of Y values.

  #store x, y, dmin for each Point
  points = zeros(n,3)

  check_dist = true
  if s == 0
    if m == 0
      check_dist = false #can put points arbitrarily close so may as well not check
      dmin = 0
    else
      d = [m] #julia can not draw from a normal distribution with sd 0
    end
  else
    d = Normal(m, s) #draw distances from normal distribution
  end

  i = 1; trials = 0
  while i <= n
    add = true
    coords = rand(2,1) #two floats [0:1]
    coords = coords.*[xhi-xlo, yhi-ylo]+[xlo, ylo] #convert to range
    #println(coords)
    if check_dist #only check distance if not m, s = 0, 0
      dmin = rand(d,1,1)[1]
      for j in 1:(i-1) #check distance to all other points
        if add #only check distance if we have not already found a point too close
          if norm(coords-points[j,1:2]) < max(dmin, points[j,3]) #too close to point
            add = false
            Print && println(coords, " and ", points[j,:], " too close with dmin ", dmin)
            trials += 1
            if trials > thresh
              #may be impossible/extremely unlikely to add new point so we stop
              println("Have failed to add ", thresh,
                      " points in a row. Exiting. Added ", i-1)
              return(zeros(0,0))
            end
          end
        end
      end
    end
    if add #add point to array
      if hs
        points[i, 1:3] = hcat(coords', dmin)
      else
        points[i, 1:3] = hcat(coords', 0)
      end
      i += 1
      trials = 0 #reset number of trials
    end
  end
  if Plot
    plot_points(points, pointsize = (m/2), filename="dmin",
              xlo = xlo, xhi = xhi, ylo = ylo, yhi = yhi)
  end
  return(points)
end
