
using Random, Distributions, PyCall, PyPlot, LinearAlgebra, DelimitedFiles
using Plots
@pyimport matplotlib.patches as patch

println("\n\nNew simulation")

function plot_points(points; pointsize=10, xlo=0, xhi=400, ylo=0, yhi=400,
                    xlab="", ylab="", title="", filename="default")

  figsize = [5,5] .* [yhi-ylo, xhi-xlo] ./ mean([yhi-ylo, xhi-xlo])
  #figure(figsize = figsize)
  fig, ax = subplots(figsize = figsize)
  xlabel(xlab)
  ylabel(ylab)
  xlim([xlo,xhi])
  ylim([ylo,yhi])
  println("plotting", size(points))
  circles = []
  for i in 1:size(points)[1]
    circles = vcat(circles,
              matplotlib[:patches][:Circle](
              (points[i,1], points[i,2]), radius = pointsize, facecolor="b"))
  end
  p = matplotlib[:collections][:PatchCollection](circles)
  ax[:add_collection](p)
  println("plotted")
  #plot(xs, ys, "bo", markersize=12.8)
  savefig(filename)
  close()
end

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
    plot_points(points, pointsize = 3,
              xlo = xlo, xhi = xhi, ylo = ylo, yhi = yhi)
  end
  return(points)
end

function calc_RI(points; Print = true)
  n = size(points)[1]

  n < 3 && return(Inf) #sd is zero in this case

  dmins = zeros(n)

  for i = 1:n
    point = points[i,1:2]
    dmin = Inf
    for j = vcat(1:i-1, i+1:n) #don't need self-distance. in julia, (1:0) is empty
      d = norm(point-points[j,1:2])
      if d < dmin
        dmin = d
      end
    end
    dmins[i] = dmin
  end
  RI = mean(dmins)/std(dmins)
  Print && println("RI is ", RI, " mean ", mean(dmins), " std ", std(dmins))
  return(RI)
end

function run_sims(;N = 1000, n=200, m=0, s=0, xlo=0, xhi=1000, ylo=0, yhi=1000)
  RIs = zeros(N)
  for i in 1:N
    points = dmin2d(n, m, s, xlo, xhi, ylo, yhi,
                    Print=false, Plot=false)
    RIs[i] = calc_RI(points, Print=false)
  end
  sort!(RIs, rev=true) #sorts in place
  println("50th RI ", RIs[50], " mean ", mean(RIs),
          " median ", median(RIs), " std ", std(RIs))
  return(RIs[50])
end

function vary_params(;ns = [20, 60, 100, 150, 200, 300, 400, 550, 700, 850, 1000],
                    ratios = [1, 3, 6, 10, 15, 20, 27.5, 35, 42.5, 50], Plot=true)
  #keep area constant, varia x:y ratio

  results = zeros(length(ns), length(ratios))

  for (i, n) in enumerate(ns)
    for (j, ratio) in enumerate(ratios)
      print("n: ", n, " ratio: ", ratio, ":    ")
      RI = run_sims(n=n, xhi = sqrt(10^6)*ratio, yhi=sqrt(10^6/ratio))
      results[i, j] = RI
    end
  end
  Plot && heatmap(results, ratios, ns, xlab="ratio", ylab="points",
            filename="RI50s.png", Title="variation of RI50")
  writedlm("vary_params.txt", results)
  return(results)
end

function heatmap(results, xs, ys; xlab="", ylab="",
        filename="default", Title="default")

  figure()
  imshow(results, cmap=ColorMap("gray_r"), vmin = minimum(results),
    vmax = maximum(results))
    print("plotted")
  xticks(0:(length(xs)-1), xs, rotation=90)
  yticks(0:(length(ys)-1), ys) #, fontstyle = 'italic', fontweight='bold')
  colorbar()
  xlabel(xlab)
  ylabel(ylab)
  title(Title)
  savefig(filename, bbox_inches = "tight")
  #show()
  close()
end


function calc_sim(m, s, ref; n = 100, xlo = 0, xhi = 400, ylo = 0, yhi = 400, Print=true)

  L = size(ref)[1]
  RIs = zeros(n); us = zeros(n)
  RIs[1] = calc_RI(ref, Print=false)
  for i in 2:n
    points = dmin2d(L, m, s, xlo, xhi, ylo, yhi,
                    Print=false, Plot=false)
    RIs[i] = calc_RI(points, Print=false)
  end
  for i in 1:n
    #u is magnitude of difference from RI to mean of other RIs
    u = abs( RIs[i] - 1/(n-1)*( sum(RIs[1:i-1])+sum(RIs[i+1:n]) ) )
    us[i] = u
  end

  Print && println("u1: ", us[1], " mean: ", mean(us), " sd: ", std(us))
  return us[1]
end

function scan_ms(;ms=0:0.5:21, ss=0:0.15:6.3, Plot=true, ref=ref, niter=1000)
  #ms = 0:10:20
  #ss = 0:0.5:5
  #niter = 10
  results = zeros(length(ms), length(ss))
  result_mat = zeros(0,3)
  for (i, m) in enumerate(ms)
    for (j, s) in enumerate(ss)
      print(m, " ", s, ": ")
      u = calc_sim(m, s, ref, n=niter)
      results[i,j]=u
      result_mat = [result_mat; [m, s, u]']
    end
  end
  writedlm("scan_results.txt", result_mat)
  writedlm("scan_result_mat.txt", results)
  Plot && heatmap(results, ss, ms, ylab="mean", xlab="standard deviation",
        filename="similarity.png", Title="variation of u-score")
  minu = findmin(results)
  println("min u is ", minu[1], " at (",
      ms[minu[2][1]], ",", ss[minu[2][2]], ")")

  return result_mat
end

function steepest_descent(; thresh=0.01, delta=0.05, rate = 1.5,
              nlim=10000, ref=ref, niter=1000)
  #write path to file so we can plot it after on top of the heatmap or as a 3d plot

  #init between 5, 15 and 1, 4

  params = rand(2).*[21, 6.3]
  print("initial params ", params, " :   ")
  u = calc_sim(params[1], params[2], ref, n=10)
  results = append!(copy(params), u)'
  err = thresh+1
  n = 0

  while u > thresh && n <= nlim
    niter = Int(round(max(5000*exp(-7*u), 20))) #adaptive n
    delta = u/2
    rate = u/1
    us = [calc_sim(params[1]+delta, params[2], ref, Print=false, n=niter),
          calc_sim(params[1], params[2]+delta, ref, Print=false, n=niter)]
    params += ([u, u] - us)/delta*rate
    print("new params ", round.(params, digits=3), " n=", niter, " :   ")
    u = calc_sim(params[1], params[2], ref, n=niter)
    results = [results; append!(copy(params), u)']
    #err = abs(u-u0)
    n += 1
    #println(err, ' ', n)
  end
  writedlm("descent_results.txt", results)
  return results
end

function test_sd(;n = 10, niter = 1000, ref=ref)
  us = zeros(10)
  for i in 1:10
    us[i] = calc_sim(12, 3, ref, n=niter)
  end
  sd = std(us)
  println("sd = ", sd) #result 0.011
  return sd
end

function plot_descent()
  scan = readdlm("scan_results.txt")
  N = Int(sqrt(size(scan)[1]))
  ms = reshape(scan[:,1], N, N)
  ss = reshape(scan[:,2], N, N)
  us = reshape(scan[:,3], N, N)

  print(size(us))

  heatmap(convert(Array, us'), 0:0.3:6.3, 0:1:21, ylab="mean", xlab="standard deviation",
        filename="similarity_gray.png", Title="variation of u-score")

  descent = readdlm("descent_results.txt")
  figure()
  plot_surface(ms, ss, us, cmap=ColorMap("gray_r") )
  p = plot(descent[:,1], descent[:,2],
      descent[:,3].+0.3, "k-", linewidth=2)

  savefig("projected_descent.png")
  #show()
  close()
end


function plot_packings()

  xs=repeat(0:20:400, outer=21)
  ys=repeat(0:20:400, inner=21)
  plot_points(hcat(xs, ys), filename = "square.png", pointsize=10,
      xlo = 0, xhi = 400, ylo = 0, yhi = 400)

  xs = vcat(repeat(0:20:400, outer=12), repeat(10:20:390, outer=12))
  ys = vcat(repeat(0:2*sqrt(300):400, inner=21),
            repeat(sqrt(300):2*sqrt(300):400, inner=20))
  plot_points(hcat(xs, ys), filename = "hexagonal.png", pointsize=10,
      xlo = 0, xhi = 400, ylo = 0, yhi = 400)

end

function saturation(n;N=10, method = "seq")
  for i in 1:N
    if method == "seq"
      trial = size(dmin2d(n, 20, 0, 0, 400, 0, 400, Plot=false))[1]
    else
      trial = size(birth_death_model(n, m=20, xlo=0, xhi=400, ylo=0, yhi=400))[1]
    end
    if trial == n
      return false #not saturated yet
    end
  end
  return true #failed N times in a row; saturated
end

function test_max(;nmax = 492, method = "seq")
  #do this imperically; if we fail 10 consecutive times we consider it failed
  nmin = 0
  while (nmax-nmin) > 1
    global n = Int(round(mean([nmin, nmax])))
    println("new n: ", n, "  nmin: ", nmin, "  nmax: ", nmax, " ")
    global sat = saturation(n, method=method)

    if sat
      nmax = n
      println("saturated")
    else
      nmin = n
      println("unsaturated")
    end
  end
  if sat
    return (n-1) #couldn't add this many points
  else
    return n #could just add this many points
  end
end

function repeat_max(N = 10; method = "seq")
  ns = zeros(N)
  for i in 1:N
    ns[i] = test_max(method = method)
  end
  println("\nn's are:", ns)
  return ns
end


function birth_death_round(points, m, xlo, xhi, ylo, yhi)
  #I don't entirely understand what we're meant to do here...
  n = size(points)[1]
  sequence = randperm(n) #numbers 1:n in random order
  nmax = 10000
  #println(points)

  for i in sequence
    ## Point i must now be killed, and a new point
    ## positioned (born) randomly subject to satisfying
    ## the minimal distance constraint.
    #YOUR CODE HERE
    niter = 0
    add = false
    while (!add) & (niter < nmax)
      #println("trying to add point")
      niter += 1
      global coords = rand(2,1).*[xhi-xlo, yhi-ylo]+[xlo, ylo]
      add = true
      for j in vcat(1:(i-1), (i+1):n) #check distance to all other points
        if add
          if norm(coords-points[j,1:2]) < m #too close to point
            add = false
          end
        end
      end
    end
    #println(coords, add)
    if add
      points[i,:] = coords'
    else
      println("Have failed to add ", nmax,
              " points in a row. Exiting. Added ", i-1)
      return(zeros(0,0))
    end
  end
  return points
end

function birth_death_model(n; m=20, xlo=0, xhi=400, ylo=0, yhi=400, nepochs = 10)
  #we hardcode here s=0
  d = [m]
  points = zeros(n,2)
  points[:,1] = (rand(n).*(xhi-xlo)).+xlo
  points[:,2] = (rand(n).*(yhi-ylo)).+ylo
  for epoch = 1:nepochs
    points = birth_death_round(points, m, xlo, xhi, ylo, yhi )
    #println(size(points))
  end
  return points
end


####Now run actual code

points = dmin2d(200, 30, 5, 200, 1000, 100, 900, Print=false, Plot=false)
#RI = calc_RI(points)
#results = vary_params()
ref = readdlm("spa3_real.dat", Float64)
#plot_packings()
#calc_sim(10, 5, ref)
#sd = test_sd()
#result_mat = scan_ms(ref=ref)
#results = steepest_descent()
plot_descent()
#n = test_max()
#a = birth_death_model(20)
#ns = repeat_max()
#ns = repeat_max(method="bd")
#writedlm("nmax_list_bd.txt", ns)
