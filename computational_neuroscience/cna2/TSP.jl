#Code for the travelling salesman problem
#Durbin R, Willshaw D (1987)
#An analogue approach to the travelling salesman problem using an elastic net method
#Nature 326:689–691.

using PyCall, PyPlot, LinearAlgebra, Distributions,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats

Random.seed!(07031106)

function swap_cities(cities, N, L0)
    cities0 = copy(cities)
    i, j = Rand(1:N, 2)

    if j > i inds = [i+1:j-1]
    elseif i > j inds - [i+1:N ; 1:j-1]
    else inds = [1] end

    cities[inds,:] = cities[reverse(inds),:]
    L = L0 +
        norm(cities[mod(i+1-1, N)+1,:]-cities[i,:]) +
        norm(cities[j,:]-cities[mod(j-1-1)+1,:]) -
        norm(cities0[mod(i+1-1, N)+1,:]-cities0[i,:]) -
        norm(cities0[j,:]-cities0[mod(j-1-1)+1,:])

    return cities, L
end
function sim_annealing(cities0, nmax = 10^7, T0=10)
    #reverse path between two points
    N = size(cities)[1]

    cities0 = cities0[shuffle(1:end), :]
    L0 = get_length(cities0)

    for i=n = 1:nmax
        T = T0/log(n) - T0/log(nmax)
        cities, L = swap_cities(cities0, N)
        if L < L0
            L0, cities0 = L, cities
        else
            p = exp(-(L-L0)/T)
            if rand()<p
                L0, cities0 = L, cities
            end
        end
    end
    print("final length is: ", L0)
    return cities0, L0
end


function plot_route(cities, points; name="figures/testtpm.png")
    figure()
    println(size(points))
    points = vcat(points, points[1, :]')
    println(size(points))
    plot(cities[:,1], cities[:,2], "ko")
    plot(points[:,1], points[:,2], "b-")
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
    N = size(cities)[1]; M = size(points)[1]

    ws = zeros(N,M)
    for i = 1:N #construct weight matrix
        xi = cities[i,:]
        for j = 2:M
            yj = points[j,:]
            ws[i,j] = exp( - (norm(yj-xi)^2) / (2*K^2) )
        end
        ws[i,:] /= sum(ws[i,:])
    end

    newpoints = zeros(M,2) #vector for storing new points
    for j = 1:M
        yj = points[j,:]
        #println(size(reshape(ws[:,j], N, 1)))
        #println(size(cities[:,:]))
        #println(size(yj))
        #println(size((cities[:,:].-yj')))
        delta = alpha * sum( (reshape(ws[:,j], N, 1) .* (cities[:,:].-yj')), dims=1 ) # first term
        #use modulus calculations to reflect circular structure in data
        delta = delta[1,:] + beta*K*(points[mod(j+1-1, M)+1,:]-2*yj+points[mod(j-1-1,M)+1,:]) #second term
        newpoints[j,:] = yj + delta
    end
    return newpoints
end

function get_init_points(cities, M)
    centroid = mean(cities, dims=1)[1,:]
    points = zeros(M,2)
    for i = 1:M
        ang = i/M * 2 * pi
        pos = centroid + (0.099+rand(1)[1]*0.002) .* [cos(ang); sin(ang)]
        points[i,:] = pos
    end
    return points
end

function get_length(points)
    points = vcat(points, points[1, :]')
    L = 0
    for i=2:size(points)[1]
        L += norm(points[i,:]-points[i-1,:])
    end
    println("length: ", L)
    return L
end

function TSP(; cities=cities, K = 0.2, alpha=0.2, beta=2.0, Klim=0.001,
            name="test_tsp", decay=0.0005, M=3, Plot=true, plot_int = false)
    N = size(cities)[1]
    M = Int(round(M*N))
    points = get_init_points(cities, M)
    Plot && plot_route(cities, points, name="figures/"*name*"_1.png")
    n = 0
    while K > Klim
        n += 1
        if plot_int
            if n%200 == 0
                plot_route(cities, points,
                    name="figures/int/"*name*"_K"*string(K)*".png")
            end
        end
        #if n%25 == 0 K *= decay end
        K -= decay*K
        points = update_round(cities, points, K, alpha, beta)
    end
    println(n, " ", K)
    #println(points)
    Plot && plot_route(cities, points, name="figures/"*name*"_2.png")
    Plot && get_length(points)
    return cities, points
end

function test_parameter_space()
    cities = readdlm("cities.dlm")
    #cities = readdlm("cities_small.dlm")
    decays = [0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001]
    Ms = [2, 3, 4, 5, 10, 50]
    decays = [0.01, 0.005, 0.001, 0.0005, 0.0001]
    Ms = [1.2, 1.5, 2,3,5]
    Ls = zeros(length(decays), length(Ms))
    Ts = zeros(length(decays), length(Ms))
    for (i, decay) in enumerate(decays)
        for (j, M) in enumerate(Ms)
            t = time()
            cities, points = TSP(cities=cities, decay = decay, M = M, name="scan", Plot=false)
            T = time()-t
            Ts[i,j] = T
            L = get_length(points)
            println("(", decay, " ", M, ") ", L, " ", T)
            Ls[i,j] = L
        end
    end
    writedlm("Ls.dlm", Ls)
    writedlm("Ts.dlm", Ts)
    writedlm("Ms.dlm", Ms)
    writedlm("decays.dlm", decays)
    heatmap(Ls, Ms, decays; xlab="M", ylab="lambda",
        filename="figures/Ls", Title="default")
    heatmap(Ts, Ms, decays; xlab="M", ylab="lambda",
        filename="figures/Ts", Title="default")
    return Ms, decays, Ls, Ts
end

function compare_methods(;n = 30, name = "rand1")
    if name == "orig"
        cities = readdlm("cities.dlm")
    else
        cities = rand(100,2)
        writedlm("cities_"*name*".dlm", cities)
    end
    t = time()
    cities, points = TSP(name="compare_elastic_"*name, decay=0.0005, M=1.5, plot_int=false)
    Torig = time()-t
    Lorig = get_length(points)

    Ls = zeros(n)
    Ts = zeros(n)
    for i=1:n
        t = time()
        cities0, L0 = sim_annealing(cities)
        Ts[i] = time()-t
        Ls[i] = L0
    end
    println("Lorig: ", Lorig, "   Torig: ", Torig)
    println("Mean: ",mean(Ls),"  std: ",std(Ls),"  min: ",minimum(Ls),"  max: ",maximum(Ls))
    println("Mean: ",mean(Ts),"  std: ",std(Ts),"  min: ",minimum(Ts),"  max: ",maximum(Ts))
end

cities_small = rand(10,2)
writedlm("cities_small.dlm", cities_small)
#cities = rand(100,2)
#writedlm("cities.dlm", cities)
cities = readdlm("cities.dlm")

#cities, points = TSP(name="test_tsp", decay=0.0005, M=1.5, plot_int=true)
#test_parameter_space()


#points = get_init_points(cities, 25)
#plot_route(cities, points, name="figures/test1.png")
#for i = 1:10000
#    global points = update_round(cities, points, 0.2, 0.2, 2.0)
#end
#plot_route(cities, points, name="figures/test2.png")
