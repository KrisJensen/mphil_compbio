#Gillespie??

using PyCall, PyPlot, LinearAlgebra, Distributions,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats

function rates_ass(x; alpha=1, lambda=1, tau1=1, tau2=2, tau3=4)
    #given a state vector x returns the rates of our 6 reactions
    #use different function here for different system
    return [alpha #x1 -> x1+1
            x[1] / tau1 #x1 -> x1-1
            lambda * x[1] #x2 -> x2+1
            x[2] / tau2 #x2 -> x2-1
            lambda * x[1] #x3 -> x3+1
            x[3] / tau3] #x3 -> x3-1
end

function gillespie(x0, deltas; niter = 1000, lambda=1, alpha=1,
                    tau1=1, tau2=2, tau3=4)
    #runs gillespi algorithm for niter iterations given parameters
    #x0 is a column vector of the initial state.
    #deltas is a list of state changes corresponding to the different reactions

    t = 0 #keep track of time
    x = x0 #colvector
    Nx = length(x) #number of variables
    Ndel = length(deltas) #number of reactions
    xs = zeros(Nx, niter+1); xs[:,1] = x0 #keep track of state
    rxns = zeros(Ndel, niter) #store which reaction occurs
    ts = zeros(niter+1) #keep track of time
    dts = zeros(niter) #store vector of delta ts (redundant)
    xweights = zeros(50000,Nx) #occupancy of each x value

    for n in 1:niter
        rates = rates_ass(x, lambda=lambda, alpha=alpha,
                        tau1=tau1, tau2=tau2, tau3=tau3) #calculate new rates
        rtot = sum(rates) #total rate
        dt = -log(rand())/rtot
        for k in 1:Nx xweights[x[k]+1, k] += dt end #store occupancy

        i = sample(1:Ndel, Weights(rates ./ rtot)) #draw reaction to occur (int in 1:6)
        t += dt #update time
        x += deltas[i] #update state
        xs[:,n+1] = x #store state
        ts[n+1] = t #store time
        dts[n] = dt #store dt
        rxns[i,n] += 1 #update reaction that occurred
    end
    xweights = xweights ./ t #normalize occupancy for probability distribution
    return dts, ts, xs, xweights, rxns
end

function plot_flux(rxns, xs, dts, ts)
    #plot how the normalized difference in flux changes over time
    #rxns is a matrix of reactions that occurred

    xs = xs[:, 1:(end-1)] .* reshape(dts, 1, length(dts)) #scale xs by time
    coords = 1000:1000:length(dts) #timepoints to plot
    N = length(coords)
    normfluxes = zeros(3, N) #array for storing fluxes
    for (n, coord) in enumerate(coords)
        #for every timepoint, find the mean normalized flux up to this timepoint
        for i in 1:3 #do this for x1, x2, x3
            t = ts[coord]
            Rip = sum(rxns[2*i-1,1:coord])/t #positive flux
            Rim = sum(rxns[2*i,1:coord])/t #negative flux
            xmean = sum(xs[i,1:coord])/t #mean of x
            normfluxes[i,n] = (Rip-Rim)/xmean #normalized flux
        end
    end
    println("\nfinal fluxes: ", normfluxes[:, end], "\n")
    figure() #plot fluxes

    cols = ["r", "g", "b"]
    for i in 1:3
        plot(coords, normfluxes[i,:], cols[i]*"--")
    end
    plot([0; coords[end]], [0; 0], "k--")
    legend(["x1", "x2", "x3"])
    savefig("figures/testfluxes.png")
    close()

end


function get_noise(dts, xs; tau1=1, tau2=2, tau3=4, Print=true)
    #return sigma_x^2 / <x>^2 for the different reactions
    N = size(xs)[1] #number of parameters
    x2s = xs.^2
    ttot = sum(dts)
    #initialize arrays for storing noise terms
    noises = zeros(N); means = zeros(N)
    extrinsic = zeros(N); intrinsic = zeros(N)
    for i in 1:3 #iterate through variables
        #mean x value over time
        meanx = sum(xs[i, 1:(end-1)] .* dts)/ttot
        #mean x^2 value over time
        meanx2 = sum(x2s[i, 1:(end-1)] .* dts)/ttot
        Print && println(meanx2," ", meanx)
        varx = meanx2 - meanx^2 #variance
        noise = varx/(meanx^2) #CV^2
        noises[i] = noise #store CV^2
        means[i] = meanx #store mean value
        intrinsic[i] = 1/meanx #intrinsic noise (poisson loss)
        if i == 1
            expec = 1/meanx #only poisson
            extrinsic[i] = 0 #no extrinsic noise
        elseif i == 2
            #calculate expected value from part A
            expec = 1/meanx+tau1/((tau1+tau2)*means[1])
            extrinsic[i] = tau1/((tau1+tau2)*means[1]) #extrinsic noise
        elseif i == 3
            #calculate expected value from part B
            expec = 1/meanx+tau1/((tau1+tau3)*means[1])
            extrinsic[i] = tau1/((tau1+tau3)*means[1]) #extrinsic noise
        else
            println("something is wrong...")
        end
        Print && println(i, "  noise: ", noise, "  exp: ", expec)
    end

    eta = zeros(3,3) #calculate covariance matrix
    for i = 1:3
        for j = i:3 #find all elements
            xixj = xs[i, :] .* xs[j, :]
            meanxi = sum(xs[i, 1:(end-1)] .* dts)/ttot #<xi>
            meanxj = sum(xs[j, 1:(end-1)] .* dts)/ttot #<xj>
            meanxixj = sum(xixj[1:(end-1)] .* dts)/ttot #<xixj>
            cov = meanxixj - meanxi*meanxj
            cov /= (meanxi * meanxj) #cov/<xi><xj>
            eta[i,j] = cov
            eta[j,i] = cov
        end
    end

    return intrinsic, extrinsic, noises, eta #return var/mean^2
end

function plot_figs(ts, xs, xweights; fname="figures/test")
    #plot part of trajectory and histogram of x values
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

function vary_params(param, values; fname="figures/test_noise",
                    tau1=3, tau2=2, tau3=4, lambda=1, alpha=10, niter=500000)
    #calculate & plot intrinsic and extrinsic noise for a range of parameters
    #param is name of parameter to vary, values is set of values
    N = length(values) #number of parameter values
    exts = zeros(3, N); ints = zeros(3, N) #store extrinsic and intrinsic noise
    errs = zeros(3,N) #store relative errors
    for (i, val) in enumerate(values)
        #find the parameter to vary
        if param == "lambda" lambda = val
        elseif param == "alpha" alpha = val
        elseif param == "tau1" tau1 = val
        else println("some is wrong...") end
        #run Gillespie simulation
        x0 = [Int(round(alpha*tau1));
            Int(round(alpha*tau1*tau2*lambda));
            Int(round(alpha*tau1*tau3*lambda))]
        dts, ts, xs, xweights, rxns = gillespie(x0, step_size, niter=niter,
                                    alpha=alpha, lambda=lambda,
                                    tau1=tau1, tau2=tau2, tau3=tau3)
        #calculate noise terms
        intrinsic, extrinsic, noises =
                    get_noise(dts, xs; tau1=tau1, tau2=tau2, tau3=tau3, Print=false)
        for j in 1:3 #store intrinsic and extrinsic noise
            ints[j, i] = intrinsic[j]
            exts[j, i] = extrinsic[j]
            #also calculate flux deviation from 0
            errs[j, i] = (sum(rxns[2*j-1,:]) - sum(rxns[2*j,:])
                        ) / sum(xs[:, 1:(end-1)] .* reshape(dts, 1, length(dts)))
        end
        println(param, " ", i, " ", val, "  err: ", errs[:,i])
    end
    figure(figsize = (5,3.5)) #plot intrinsic and extrinsic noise
    cols = ["r", "b", "g"]
    for i in 2:3
        plot(values, ints[i,:], cols[i]*"-")
        plot(values, exts[i,:], cols[i]*"--")
        plot(values, ints[i,:]./exts[i,:], cols[i]*":")
    end
    legend(["x2 intrinsic", "x2 extrinsic", "x2 int/ext",
            "x3 intrinsic", "x3 extrinsic", "x3 int/ext"])
    xscale("log"); yscale("log")
    xlabel(param); ylabel("variance")
    savefig(fname*"_noise.png", bbox_inches="tight", dpi=120)
    close()

    xmin = minimum(errs); xmax=maximum(errs)
    #bins = xmin : (xmax-xmin)/199 : xmax
    for i in 1:3
        figure()
        PyPlot.plt[:hist](errs[i,:], 20, color=cols[i])
        savefig(fname*"_hist"*string(i)*".png", bbox_inches="tight", dpi=120)
        close()
    end

end
println("\nnew simulation")


step_size = [ [1; 0; 0], [-1; 0; 0], #x1 +- 1
            [0; 1; 0], [0; -1; 0], #x2 +- 1
            [0; 0; 1], [0; 0; -1]] #x3 +- 1

tau2=2 #fixed
tau3=4 #fixed

tau1=3
lambda=1
alpha=10
niter = 100000

x0 = [0; 0; 0] #initial conditions
x0 = [alpha*tau1; alpha*tau1*tau2*lambda; alpha*tau1*tau3*lambda]


dts, ts, xs, xweights, rxns = gillespie(x0, step_size, niter=niter,
                            alpha=alpha, lambda=lambda,
                            tau1=tau1, tau2=tau2, tau3=tau3)

println(get_noise(dts, xs, tau1=tau1, tau2=tau2, tau3=tau3))
plot_figs(ts, xs, xweights)
plot_flux(rxns, xs, dts, ts)

#find noise vs. lambda
vary_params("lambda", 10 .^ (-2:(4/99):2), fname="figures/test_noise" )

#find noise vs. lambda
#vary_params("lambda", 10 .^ (-2:(4/99):2), fname="figures/lambda_noise.png" )

#find noise vs. alpha
#vary_params("alpha", 10 .^ (-2:(3/99):1), fname="figures/alpha_noise.png" )

#find noise vs. tau1
#vary_params("tau1", 10 .^ (-2:(5/99):3), fname="figures/tau1_noise.png" )
