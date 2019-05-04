#Gillespie??

using PyCall, PyPlot, LinearAlgebra, Distributions,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats

function rates_pois(x)
    #x = x, y
    lambda1 = 10
    beta1 = 1
    lambda2 = 5
    beta2 = 2
    lambda3 = 3
    beta3 = 3
    return [lambda1
            beta1*x[1]
            lambda2
            beta2*x[2]
            lambda3
            beta3*x[3]]
end

function rates_ass(x; alpha=1, lambda=1, tau1=1, tau2=2, tau3=4)

    return [alpha #x1 -> x1+1
            x[1] / tau1 #x1 -> x1-1
            lambda * x[1] #x2 -> x2+1
            x[2] / tau2 #x2 -> x2-1
            lambda * x[1] #x3 -> x3+1
            x[3] / tau3] #x3 -> x3-1
end





function gillespie(x0, deltas; niter = 1000, lambda=1, alpha=1,
                    tau1=1, tau2=2, tau3=4)
    #x0, deltas are colvectors
    t = 0
    x = x0 #colvector
    Nx = length(x)
    Ndel = length(deltas)
    #xs = zeros(length(x), 1); xs[:,1] = x0
    xs = zeros(Nx, niter+1); xs[:,1] = x0
    ts = zeros(niter+1)
    dts = zeros(niter)
    dt = 0
    xcounts = zeros(1000,Nx)
    xweights = zeros(1000,Nx)
    for k in 1:Nx xcounts[x[k]+1, k] += 1 end
    #ts = zeros(1, 1); ts[1,1] = t
    #while (t < tmax)
    for n in 1:niter
        for k in 1:Nx xweights[x[k]+1, k] += dt end #previous state
        #rates = rates_pois(x) #calculate new rates
        rates = rates_ass(x, lambda=lambda, alpha=alpha,
                        tau1=tau1, tau2=tau2, tau3=tau3) #calculate new rates
        rtot = sum(rates) #total rate
        dt = -log(rand())/rtot
        i = sample(1:Ndel, Weights(rates ./ rtot)) #draw reaction to occur (int in 1:6)

        t += dt #update time
        x += deltas[i] #update state
        for k in 1:Nx xcounts[x[k]+1, k] += 1 end
        #ts = [ts t] #store new time
        #xs = [xs x] #store new state
        xs[:,n+1] = x
        ts[n+1] = t
        dts[n] = dt
    end
    xweights = xweights ./ t
    return dts, ts, xs, xcounts, xweights
end

function get_noise(dts, xs; tau1=1, tau2=2, tau3=4, Print=true)
    #return sigma_x / <x>
    N = size(xs)[1]
    x2s = xs.^2
    ttot = sum(dts)
    noises = zeros(N)
    means = zeros(N)
    extrinsic = zeros(N); intrinsic = zeros(N)
    for i in 1:3
        meanx = sum(xs[i, 1:(end-1)] .* reshape(dts, length(dts), 1))/ttot
        meanx2 = sum(x2s[i, 1:(end-1)] .* reshape(dts, length(dts), 1))/ttot
        Print && println(meanx2," ", meanx)
        varx = meanx2 - meanx^2
        normstd = sqrt(varx)/meanx
        #normstd = std(xs[i,:])/meanx
        noises[i] = normstd
        means[i] = meanx
        intrinsic[i] = 1/meanx
        if i == 1
            expec = 1/meanx
            extrinsic[i] = 0
        elseif i == 2
            expec = 1/meanx+tau1/((tau1+tau2)*means[1])
            extrinsic[i] = tau1/((tau1+tau2)*means[1])
        elseif i == 3
            expec = 1/meanx+tau1/((tau1+tau3)*means[1])
            extrinsic[i] = tau1/((tau1+tau3)*means[1])
        else
            println("something is wrong...")
        end
        Print && println(i, "  normstd: ", normstd, "  exp: ", sqrt(expec))
    end
    return intrinsic, extrinsic, noises #var/mean^2 for intrinsic/extrinsic, std/mean for total
end

function plot_figs(ts, xs, xweights; fname="figures/test")

    figure() #plot part of the trajectory
    t1, t2 = 1000, 1500
    plot(ts[t1:t2], xs[1,t1:t2], "r--")
    plot(ts[t1:t2], xs[2,t1:t2], "g--")
    plot(ts[t1:t2], xs[3,t1:t2], "b--")
    savefig(fname*"_trajec.png")
    close()

    xvals = 0:20
    figure() #plot a weighted histogram (xweights)
    cols = ["r", "g", "b"]
    for i = 1:3
        lambda = lambdas[i]; beta = betas[i]
        plot(0:(size(xweights)[1]-1), xweights[:,i], cols[i]*"-")
    end
    legend(["x1", "x2", "x3"])
    xlim(0, 20)
    savefig(fname*"_hist.png")
    close()
end

function vary_params(param, values; fname="figures/test_noise.png",
                    tau1=1, tau2=2, tau3=4, lambda=1, alpha=1, niter=100000)
    #param is name of parameter to vary, values is set of values
    N = length(values)
    exts = zeros(3, N); ints = zeros(3, N)
    for (i, val) in enumerate(values)
        println(param, " ", i, " ", val)
        if param == "lambda" lambda = val
        elseif param == "alpha" alpha = val
        elseif param == "tau1" tau1 = val
        else println("some is wrong...") end
        dts, ts, xs, xcounts, xweights = gillespie(x0, step_size, niter=niter,
                                    alpha=alpha, lambda=lambda,
                                    tau1=tau1, tau2=tau2, tau3=tau3)
        intrinsic, extrinsic, noises =
                    get_noise(dts, xs; tau1=tau1, tau2=tau2, tau3=tau3, Print=false)
        for j in 1:3
            ints[j, i] = intrinsic[j]
            exts[j, i] = extrinsic[j]
        end
    end
    figure(figsize = (5,3.5))
    plot(values, ints[2,:], "b-")
    plot(values, exts[2,:], "b--")
    plot(values, ints[3,:], "g-")
    plot(values, exts[3,:], "g--")
    legend(["x2 intrinsic", "x2 extrinsic", "x3 intrinsic", "x3 extrinsic"])
    xscale("log")
    yscale("log")
    xlabel("lambda")
    ylabel("variance")
    savefig(fname, bbox_inches="tight", dpi=120)
    close()
end
println("\nnew simulation")

x0 = [0; 0; 0] #initial conditions
step_size = [ [1; 0; 0], [-1; 0; 0], #x1 +- 1
            [0; 1; 0], [0; -1; 0], #x2 +- 1
            [0; 0; 1], [0; 0; -1]] #x3 +- 1

tau2=2 #fixed
tau3=4 #fixed

tau1=1
lambda=1
alpha=1
niter = 100000


dts, ts, xs, xcounts, xweights = gillespie(x0, step_size, niter=niter,
                            alpha=alpha, lambda=lambda,
                            tau1=tau1, tau2=tau2, tau3=tau3)

println(get_noise(dts, xs, tau1=tau1, tau2=tau2, tau3=tau3))

plot_figs(ts, xs, xweights)

vary_params("lambda", 10 .^ (-2:0.1:1) )
