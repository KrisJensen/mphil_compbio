#code for using the elastic net algorithm to model ocular ocular dominance

using PyCall, PyPlot, LinearAlgebra, Distributions,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats

Random.seed!(21031633) #random seed for reproducibility

function plot_pref(cities, points; name = "figures/occ2d")
    #function for plotting results of a simulation
    N = size(cities)[2]
    M = size(points)[1]

    #plot ocular dominance
    zval = points[:, :, 3]
    zval = sign.(zval)
    figure()
    imshow(zval, cmap=ColorMap("gray_r"), vmin = -1,
            vmax = 1)
    xticks([], [])
    yticks([], [])
    savefig(name*"_dominance.png", bbox_inches = "tight")
    close()

    #plot x mapping
    xval = points[:, :, 1]
    figure()
    imshow(xval, cmap=ColorMap("gray_r"), vmin = minimum(xval),
            vmax = maximum(xval))
    xticks([], [])
    yticks([], [])
    savefig(name*"_xvals.png", bbox_inches = "tight")
    close()

    #plot y mapping
    yval = points[:, :, 2]
    figure()
    imshow(yval, cmap=ColorMap("gray_r"), vmin = minimum(yval),
            vmax = maximum(yval))
    xticks([], [])
    yticks([], [])
    savefig(name*"_yvals.png", bbox_inches = "tight")
    close()

    #plot all points in 3D
    figure()
    plot3D(
        reshape(cities[1,:,:,1], N*N),
        reshape(cities[1,:,:,2], N*N),
        reshape(cities[1,:,:,3], N*N), color=[0,0,1,0.5],
        linestyle = "", marker="o", markersize = 1)
    plot3D(
        reshape(cities[2,:,:,1], N*N),
        reshape(cities[2,:,:,2], N*N),
        reshape(cities[2,:,:,3], N*N), color=[1,0,0,0.5],
        linestyle = "", marker="o", markersize = 1)
    plot3D(
        reshape(points[:,:,1], M*M),
        reshape(points[:,:,2], M*M),
        reshape(points[:,:,3], M*M), color=[1,0,1,0.5],
        linestyle = "", marker="o", markersize = 1)
    savefig(name*"_positions.png", bbox_inches = "tight")
    close()
end

function update_round(cities, points, K, alpha, beta, Print=false)
    #performs one update step of the elastic net
    N = size(cities)[2]; M = size(points)[1]
    ws = zeros(2,N,N,M,M)

    #construct weight matrix
    for k = 1:2 #right/left eye
        for a = 1:N #retinal position
            for b = 1:N
                xab = cities[k,a,b,:]
                for i = 1:M  #cortical position
                    for j = 1:M
                        yij = points[i,j,:]
                        ws[k,a,b,i,j] = exp( - (norm(yij-xab)^2) / (2*K^2) )
                        if isnan(ws[k,a,b,i,j])
                            return 0, 0 #things are broken
                        end
                    end
                end
                Print && println("\n")
                Print && println(ws[k,a,b,:,:])
                ws[k, a, b, :, :] /= sum(ws[k,a,b,:,:]) #normalize
                Print && println(ws[k,a,b,:,:], "    ", sum(ws[k,a,b,:,:]))
            end
        end
    end

    newpoints = zeros(M,M,3) #vector for storing new points
    for i = 1:M
        for j = 1:M
            yij = points[i, j,:] #old position
            delta = zeros(3) #change in position
            for k = 1:2
                for a = 1:N
                    for b = 1:N
                        xab = cities[k,a,b,:] #retinal point
                        delta = delta + alpha * ws[k,a,b,i,j] * (xab - yij)
                    end
                end
            end
            Print && println("\n")
            Print && println(delta)

            #add elastic contributions
            n = 0
            if (j < M) delta = delta + beta*K*points[i, j+1,:]; n+=1 end
            if (i < M) delta = delta + beta*K*points[i+1, j,:]; n+=1 end
            if (j > 1) delta = delta + beta*K*points[i, j-1,:]; n+=1 end
            if (i > 1) delta = delta + beta*K*points[i-1, j,:]; n+=1 end
            delta = delta - n*beta*K*yij

            Print && println(delta)
            newpoints[i, j,:] = yij + delta
        end
    end
    return newpoints
end



function run_sim(cities, cortex, N, M; K = 0.1, alpha=0.2, beta=2.0, Klim=0.001,
            name="occ2d_test", decay=0.001, Plot=true, plot_int = false)
    #performs an elastic net optimization of the TSP decreasing K exponentially from K
    #with decay rate decay. M specifies the number of points on the path
    #as a multiple of N

    Plot && plot_pref(cities, cortex, name="figures/"*name*"_1")
    n = 0 #number of iterations
    while K > Klim
        n += 1
        if n % 10 == 0
            println("n: ", n, "   K: ", K, "  Klim:", Klim)
            if isnan(minimum(cortex)) #simulation diverged
                println("nan encountered, exiting")
                return cities, cortex
            end
            if n % 100 == 0 #plot intermediate points
                if plot_int
                    plot_pref(cities, cortex, name="figures/int/"*name*"_K"*string(K))
                end
            end
        end
        K -= decay*K #exponential decay
        cortex = update_round(cities, cortex, K, alpha, beta) #update path
    end
    println(n, " ", K)
    Plot && plot_pref(cities, cortex, name="figures/"*name*"_2")
    return cities, cortex
end

function get_retinae(;N=20, M=40, l =0.075, d = 0.05/2)
    #N^2 is number of points in retinal planes
    #M^2 in cortex plane
    #l is inter-plane separation
    #b is inter-points separation

    #first index specifies retina, next two xy index of point, then coordinates
    cities = zeros(2,N,N,3)
    cortex = zeros(M, M, 3)
    for k = 1:2 #right vs left retina
        for i = 1:N #x location
            for j = 1:N #y location
                cities[k,i,j,:] = [i*d; j*d; l*(3-2*k)]
            end
        end
    end

    #initialize cortical positions uniformly at random in xy
    #uniform distribution in z plane
    for i = 1:M
        for j = 1:M
            cortex[i,j,:] = [rand()*N*d; rand()*N*d; rand()*0.01-0.01/2]
        end
    end

    return cities, cortex
end

#specify parameters (determined empirically)
alpha = 0.2
beta = 2
K = 0.12
N = 20
M = 40
#vary ls for figures
l = 0.15/2
l = 0.2/2
l = 0.1/2
d = 0.05

cities, cortex = get_retinae(N=N, M=M, l = l, d = d) #initialize system
cities, cortex = run_sim(cities, cortex, N, M, K=K, Klim = 0.01,
                    plot_int=true, decay = 0.0005,
                    name = "N20M40l05d05") #run simulation
