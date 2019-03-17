#Code for unsupervised Hebbian learning in
#2-input-1-output systems


using PyCall, PyPlot, LinearAlgebra, Distributions,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats

Random.seed!(06032039)

function norm2d(;n=500, m=0, sd=1, rho=-0.7)
    #same mean and sd for both distributions.
    #rho defines covariance (slope).
    #return data for a 2D Gaussian distribution
    d1 = Normal(m, sd)
    x = rand(d1, n)
    d2 = Normal(m-rho*m, sqrt(1-rho^2))
    y = rand(d2, n)+rho*x
    plot(x,y, "bo")
    show()
    savefig("testfig.png")
    close()
    println(mean(x)," ", mean(y))
    println(cov(x, y))
    return [x';y']
end

function update_weights(ws, data, dt; alpha=1)
    C = cov(dat, dims=2)



    #for i = 1:size(data)[2]
    #    u = data[:,i]
    #    v = ws' * u
    #    ws = ws + dt*( v*u-alpha*v^2*ws )
    #end
end
function update_additive(w,C,n,dt)
    m = mean(w)
    finished = false
    neww = w
    neww = w + dt * (C * w - ((w'*C*n)/sum(n))*n)
    for i in 1:length(w)
        if neww[i] < 0
            #println("zeroing weight ", i)
            neww[i] = 0 #set weight to zero
            n[i] = 0 #don't alter this weight anymore
            neww = neww*m/mean(neww) #ensure sum of weights is constant
        end
    end
    println(neww)
    return neww, n
end

function simulate_system(input;tstop=10, mult = true, corr=true, dt=0.01, w=[0.001,-0.002], alpha=1, thresh = 10^(-6))
    C = zeros(2,2); n = ones(2)
    if !corr #subtract mean of input
        input = input .- mean(input,dims=2)
    end
    for i=1:2
        for j=1:2 #fill in correlation/covariance matrix
            C[i,j] = (input[i,:]' * input[j,:])/size(input)[2]
        end
    end
    println(C)
    N = Int(tstop/dt)+1
    ws = w
    w0 = w
    err = thresh+1
    while err > thresh
        if mult #multiplicative normalization
            w = w + dt * (C * w - alpha*(w'*C*w)*w)
        else
            w, n = update_additive(w,C,n,dt)
        end
        ws = hcat(ws, w)
        err = norm(w-w0)
        w0 = w
    end
    println("final vec: ", ws[:,end])
    M = fit(PCA, input)
    println("PCs: ", M.proj[:,1], " ", M.proj[:,2])
    return ws
end

function plot_dat_weights(dat, weights, name; m=0)

    figure(figsize = (4,4))
    plot(dat[1,:], dat[2,:], "bo", MarkerSize = 1)
    #plot(weights[1,:], weights[2,:], "r-")
    rat = weights[2,end]/weights[1,end]
    if rat < 0
        ys = [-10;10]*rat.+2*m
    else
        ys = [-10;10]*rat
    end
    plot([-10;10], ys, "r:")
    #plot direction of final weight vector
    xlabel(L"u_1,w_1")
    ylabel(L"u_2,w_2")
    xlim(m-3, m+3)
    ylim(m-3, m+3)
    savefig("figures/"*name*".png")
    close()

    norms = [] #Plot norms vs time
    for i=1:size(weights)[2]
        norms = [norms;norm(weights[:,i])]
    end
    plot(norms, "b-")
    xlabel("iteration")
    ylabel("|w|")
    xlim(0, length(norms))
    ylim(0,1)
    savefig("figures/"*name*"_vec.png")
    close()
end

function make_weight_fig(;basename="mul", rho=-0.7, mult=true)
    dat1 = norm2d(n=1000, m=0, rho=rho)
    dat2 = norm2d(n=1000, m=3, rho=rho)
    dat3 = norm2d(n=1000, m=3, rho=rho)
    corrs = [true, true, false]
    ms = [0;3;3]
    #set initial weights
    if mult weights = [0.001, 0.001] else weights=[0.4, 0.6] end
    for (index, dat) in enumerate([dat1, dat2, dat3])
        ws = simulate_system(dat, corr=corrs[index], mult=mult, w=weights)
        plot_dat_weights(dat, ws, basename*"_2d_sim"*string(index), m=ms[index])
    end
end


dat = norm2d(n=1000, m=0, rho=-0.7)
ws = simulate_system(dat, mult=false, w=[0.1,0.8], corr=false)
make_weight_fig(rho=-0.7, basename="test")
