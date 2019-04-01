#Code for unsupervised Hebbian learning tasks

using PyCall, PyPlot, LinearAlgebra, Distributions,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats
Random.seed!(06032039) #random seed for reproducibility

function norm2d(;n=500, m=0, sd=1, rho=-0.7, slope=1)
    #same mean and sd for both distributions.
    #rho defines covariance, slope defines slope.
    #return data for a 2D Gaussian distribution
    d1 = Normal(m, sd)
    x = rand(d1, n)
    d2 = Normal(m-rho*m, sqrt(1-rho^2))
    y = rand(d2, n)+rho*slope*x
    #plot input data
    figure()
    plot(x,y, "bo")
    xlim(min(minimum(x), minimum(y)), max(maximum(x), maximum(y)))
    ylim(min(minimum(x), minimum(y)), max(maximum(x), maximum(y)))
    savefig("testfig.png")
    show()
    close()
    println(mean(x)," ", mean(y))
    println(cov(x, y))
    return [x';y']
end


function update_additive(w,C,n,dt)
    #performs an update step for hebbian learning with subtractive normalization
    m = mean(w)
    finished = false
    neww = w
    neww = w + dt * (C * w .* n - ((w'*C*n)/sum(n))*n)
    for i in 1:length(w)
        if neww[i] < 0
            #println("zeroing weight ", i)
            neww[i] = 0 #set weight to zero
            n[i] = 0 #don't alter this weight anymore
            neww = neww*m/mean(neww) #ensure sum of weights is constant
        end
    end
    return neww, n
end

function simulate_system(input;tstop=10, mult = true, corr=true, dt=0.01, Print=true,
                w=[0.001,-0.002], alpha=1, thresh = 10^(-6), stat=false)
    #performs a full simulation using either multiplicative (mult=true)
    #or subtractive (mult=false) normalization
    #corr=true for correlation-based training, false for covariance-based
    C = zeros(2,2); n = ones(2)
    if !corr #subtract mean of input if we use a covariance based learning rule
        input = input .- mean(input,dims=2)
    end
    for i=1:2
        for j=1:2 #fill in correlation/covariance matrix
            C[i,j] = (input[i,:]' * input[j,:])/size(input)[2]
        end
    end
    Print && println(C)
    M = fit(PCA, input) #perform PCA
    #stat=true to initialize at stationary point
    if stat w = [1, (C[1,1]-C[1,2])/(C[2,2]-C[1,2])]; w = w ./ sum(w) end
    Print && println(w)
    N = Int(tstop/dt)+1

    ws = w
    w0 = w
    err = thresh+1

    while err > thresh #thresh is weight change for convergence
        if mult #multiplicative normalization
            w = w + dt * (C * w - alpha*(w'*C*w)*w)
        else
            w, n = update_additive(w,C,n,dt) #additive
        end
        ws = hcat(ws, w) #store weights
        err = norm(w-w0)
        w0 = w
    end
    Print && println("final vec: ", ws[:,end])
    Print && println("PCs: ", M.proj[:,1], " ", M.proj[:,2])
    return ws
end

function plot_dat_weights(dat, weights, name; m=0)
    #function for plotting data and weights. m is mean of data. assume equal means

    figure(figsize = (4,4))
    plot(dat[1,:], dat[2,:], "bo", MarkerSize = 1) #plot data
    rat = weights[2,end]/weights[1,end] #slope of weight vector
    if rat < 0 #shift line so it runs through data
        ys = [-10;10]*rat.+2*m
    else
        ys = [-10;10]*rat
    end
    plot([-10;10], ys, "r:") #plot direction of final weight vector
    xlabel(L"u_1,w_1")
    ylabel(L"u_2,w_2")
    xlim(m-3, m+3)
    ylim(m-3, m+3)
    savefig("figures/"*name*".png")
    close()

    norms = [] #Plot weight vector norms vs time
    for i=1:size(weights)[2]
        norms = [norms;norm(weights[:,i])]
    end
    figure(figsize = (4,3))
    plot(norms, "b-")
    xlabel("iteration")
    ylabel("|w|")
    xlim(0, length(norms))
    ylim(0,1)
    savefig("figures/"*name*"_vec.png", bbox_inches="tight")
    close()
end

function make_weight_fig(;basename="mul", rho=-0.7, mult=true)
    #function for reproducing TN figure
    dat1 = norm2d(n=1000, m=0, rho=rho) #generate data
    dat2 = norm2d(n=1000, m=3, rho=rho)
    dat3 = norm2d(n=1000, m=3, rho=rho)
    corrs = [true, true, false]
    ms = [0;3;3]
    #set initial weights
    if mult weights = [0.001, 0.001] else weights=[0.4, 0.6] end
    for (index, dat) in enumerate([dat1, dat2, dat3])
        #run simulation and plot result for each dataset
        ws = simulate_system(dat, corr=corrs[index], mult=mult, w=weights)
        plot_dat_weights(dat, ws, basename*"_2d_sim"*string(index), m=ms[index])
    end
end

function repeat_sub(;n = 100, name="test", rho=0.7, slope=1)
    #repeatedly run subtractive normalization for data with
    #orrelation rho and slope slope. Plot histogram of results
    thetas = [] #angles from 1st PC
    for i = 1:n
        dat = norm2d(n=1000, m=0, rho=rho, slope=slope) #generate data
        winit = rand(2); winit /= sum(winit) #random initial vectors
        ws = simulate_system(dat, mult=false, w=winit, corr=false,
            stat=false, Print=false)
        w = ws[:,end] #get final vector
        M = fit(PCA, dat)
        pc1 = M.proj[:,1]
        theta = min(acos(pc1' * w / (norm(pc1)*norm(w))),
                    acos(-pc1' * w / (norm(pc1)*norm(w)))) #compare to 1st PC
        thetas = [thetas; theta]
    end
    #plot histogram of results
    figure(figsize = (4,3))
    PyPlot.plt[:hist](thetas)
    xlabel("theta")
    ylabel("frequency")
    savefig("figures/"*name*"_hist.png", bbox_inches="tight")
    close()
end

make_weight_fig(rho=-0.7, basename="mult", mult = true) #reproduce TN figure
make_weight_fig(rho=-0.7, basename="sub", mult = false) #try subtractive norm

dat = norm2d(n=1000, m=0, rho=0.7, slope=1) #initialize at stationary point
ws = simulate_system(dat, mult=false, w=[0.5,0.5], corr=false, stat=true)
plot_dat_weights(dat, ws, "stationary"; m=0)

dat = norm2d(n=1000, m=0, rho=0.995, slope=0) #1st PC is [1,0]
ws = simulate_system(dat, mult=false, w=[0.5,0.5], corr=false, stat=false)
plot_dat_weights(dat, ws, "1_0"; m=0)

repeat_sub(;n = 100, name="slope_1", rho=0.7, slope=1) #histogram for slope of 1
repeat_sub(;n = 100, name="slope_02", rho=0.7, slope=0.2) #histogram for slope 0.2
