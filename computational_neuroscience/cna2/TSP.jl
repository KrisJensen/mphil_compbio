#Code for the travelling salesman problem
#Durbin R, Willshaw D (1987)

using PyCall, PyPlot, LinearAlgebra, Distributions,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats
Random.seed!(07031106) #random seed for reproducibility

function swap_cities(cities, N, L0)
    #this function swaps the part of a path between two random cities as
    #required by the simulated annealing algorithm
    cities0 = copy(cities)
    i, j = Rand(1:N, 2) #pick cities to swap
    #get indices of section to reverse
    if j > i inds = [i+1:j-1]
    elseif i > j inds - [i+1:N ; 1:j-1]
    else inds = [1] end

    cities[inds,:] = cities[reverse(inds),:] #reverse path
    L = L0 + #calculate new path length
        norm(cities[mod(i+1-1, N)+1,:]-cities[i,:]) +
        norm(cities[j,:]-cities[mod(j-1-1)+1,:]) -
        norm(cities0[mod(i+1-1, N)+1,:]-cities0[i,:]) -
        norm(cities0[j,:]-cities0[mod(j-1-1)+1,:])

    return cities, L
end
function sim_annealing(cities0, nmax = 10^7, T0=10)
    #runs a simulated annealing optimzation of the travelling salesman problem
    #nmax is number of iterations, T0 is initial temperature
    N = size(cities)[1]
    cities0 = cities0[shuffle(1:end), :] #random initial path
    L0 = get_length(cities0) #get initial length

    for i=n = 1:nmax
        T = T0/log(n) - T0/log(nmax) #new temperature
        cities, L = swap_cities(cities0, N) #reverse path between two cities
        if L < L0 #always accept if better
            L0, cities0 = L, cities
        else
            p = exp(-(L-L0)/T) #probability of swapping
            if rand()<p
                L0, cities0 = L, cities
            end
        end
    end
    print("final length is: ", L0) #current length
    return cities0, L0
end


function plot_route(cities, points; name="figures/testtpm.png")
    #function for plotting a set of cities and a route through them
    figure()
    points = vcat(points, points[1, :]') #circular route
    plot(cities[:,1], cities[:,2], "ko") #plot cities
    plot(points[:,1], points[:,2], "b-") #plot route
    savefig(name)
    close()
end

function heatmap(results, xs, ys; xlab="", ylab="",
    filename="default", Title="default")
    #given a matrix of z values and lists of x,y values
    #plots a heatmap
    figure()
    imshow(results, cmap=ColorMap("gray_r"), vmin = minimum(results),
            vmax = maximum(results))
    print("plotted")
    xticks(0:(length(xs)-1), xs, rotation=0)
    yticks(0:(length(ys)-1), ys)
    colorbar()
    xlabel(xlab)
    ylabel(ylab)
    title(Title)
    savefig(filename, bbox_inches = "tight")
    close()
end

function update_round(cities, points, K, alpha, beta)
    #performs one update step of the elastic net
    N = size(cities)[1]; M = size(points)[1]
    ws = zeros(N,M)
    for i = 1:N #construct weight matrix
        xi = cities[i,:]
        for j = 2:M
            yj = points[j,:]
            ws[i,j] = exp( - (norm(yj-xi)^2) / (2*K^2) ) #Gaussian attraction
        end
        ws[i,:] /= sum(ws[i,:])
    end

    newpoints = zeros(M,2) #vector for storing new points
    for j = 1:M
        yj = points[j,:] #old position
        delta = alpha * sum( (reshape(ws[:,j], N, 1) .*
                (cities[:,:].-yj')), dims=1 ) #first update term
        #use modulus calculations to reflect circular structure in data
        delta = delta[1,:] + beta*K*(points[mod(j+1-1, M)+1,:]-2*yj +
                points[mod(j-1-1,M)+1,:]) #second update term
        newpoints[j,:] = yj + delta #update point
    end
    return newpoints
end

function get_init_points(cities, M)
    #get a set of M initial points spread evenly around the centroid of N cities
    centroid = mean(cities, dims=1)[1,:]
    points = zeros(M,2) #initialize array
    for i = 1:M
        ang = i/M * 2 * pi #evenly distributed angles
        pos = centroid + (0.099+rand(1)[1]*0.002) .* [cos(ang); sin(ang)] #+noise
        points[i,:] = pos
    end
    return points
end

function get_length(points)
    #get the length of a TSP path specified by points
    points = vcat(points, points[1, :]') #consider circular path
    L = 0
    for i=2:size(points)[1]
        L += norm(points[i,:]-points[i-1,:]) #add length of each line segment
    end
    println("length: ", L)
    return L
end

function TSP(; cities=cities, K = 0.2, alpha=0.2, beta=2.0, Klim=0.001,
            name="test_tsp", decay=0.0005, M=3, Plot=true, plot_int = false)
    #performs an elastic net optimization of the TSP decreasing K exponentially from K
    #with decay rate decay. M specifies the number of points on the path
    #as a multiple of N
    N = size(cities)[1]
    M = Int(round(M*N)) #number of points on path
    points = get_init_points(cities, M)
    Plot && plot_route(cities, points, name="figures/"*name*"_1.png")
    n = 0 #number of iterations
    while K > Klim
        n += 1
        if plot_int
            if n%200 == 0 #plot intermediate routes
                plot_route(cities, points,
                    name="figures/int/"*name*"_K"*string(K)*".png")
            end
        end
        K -= decay*K #exponential decay
        points = update_round(cities, points, K, alpha, beta) #update path
    end
    println(n, " ", K)
    #println(points)
    Plot && plot_route(cities, points, name="figures/"*name*"_2.png")
    Plot && get_length(points)
    return cities, points
end

function test_parameter_space()
    #investigates the perform of the TSP algorithm as a function of
    #the number of points M on the path and the decay rate lambda
    cities = readdlm("cities.dlm") #read a file of stored cities
    decays = [0.01, 0.005, 0.001, 0.0005, 0.0001]
    Ms = [1.2, 1.5, 2,3,5]
    Ls = zeros(length(decays), length(Ms)) #store path lengths
    Ts = zeros(length(decays), length(Ms)) #store completion times
    for (i, decay) in enumerate(decays)
        for (j, M) in enumerate(Ms)
            t = time()
            cities, points = TSP(cities=cities, decay = decay, M = M, name="scan", Plot=false)
            T = time()-t #time for simulation
            Ts[i,j] = T #time
            L = get_length(points) #length
            println("(", decay, " ", M, ") ", L, " ", T)
            Ls[i,j] = L
        end
    end
    writedlm("Ls.dlm", Ls) #write data
    writedlm("Ts.dlm", Ts)
    writedlm("Ms.dlm", Ms)
    writedlm("decays.dlm", decays)
    heatmap(Ls, Ms, decays; xlab="M", ylab="lambda", #plot results
        filename="figures/Ls", Title="default")
    heatmap(Ts, Ms, decays; xlab="M", ylab="lambda",
        filename="figures/Ts", Title="default")
    return Ms, decays, Ls, Ts
end

function compare_methods(;n = 30, name = "rand1")
    #compare elastic net to simulated annealing for solving the TSP
    #n specifies number of simulated annealings to run
    if name == "orig"
        cities = readdlm("cities.dlm") #use our pre-defined points
    else
        cities = rand(100,2) #generate new points and save them for future use
        writedlm("cities_"*name*".dlm", cities)
    end
    t = time()
    cities, points = TSP(name="compare_elastic_"*name, decay=0.0005, M=1.5, plot_int=false)
    Torig = time()-t #time for elastic net (deterministic so only run one simulation)
    Lorig = get_length(points) #length for elastic net

    Ls = zeros(n)
    Ts = zeros(n)
    for i=1:n #try n simulations of simulated annealing
        t = time()
        cities0, L0 = sim_annealing(cities) #run simulated annealing
        Ts[i] = time()-t #store time take
        Ls[i] = L0 #store length
    end
    #print summary data
    println("Lorig: ", Lorig, "   Torig: ", Torig)
    println("Mean: ",mean(Ls),"  std: ",std(Ls),"  min: ",minimum(Ls),"  max: ",maximum(Ls))
    println("Mean: ",mean(Ts),"  std: ",std(Ts),"  min: ",minimum(Ts),"  max: ",maximum(Ts))
end

cities = rand(100,2) #generate cities
#writedlm("cities.dlm", cities) #write to file
cities = readdlm("cities.dlm") #use previously generated cities for ease of comparison

test_parameter_space() #find optimum parameters for elastic net simulation
#run single elastic net simulation
cities, points = TSP(name="test_tsp", decay=0.0005, M=1.5, plot_int=true)
compare_methods() compare elastic net and simulated annelaing
