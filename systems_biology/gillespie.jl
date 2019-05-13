#Script for implementing the Gillespie algorithm
#varying parameter values
#and plotting various diagrams used in the assignment

using PyCall, PyPlot, LinearAlgebra, Distributions, Random,
    DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats

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

function gillespie(x0, deltas; niter = 10^7, lambda=1, alpha=10,
                    tau1=3, tau2=2, tau3=4)
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

        i = sample(1:Ndel, Weights(rates ./ rtot)) #draw reaction (int in 1:6)
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

function plot_flux(rxns, xs, dts, ts; fname="figures/testfluxes.png",
                    ylimvals=[-0.02 0.02])
    #plot how the normalized difference in flux changes over time
    #rxns is a matrix of reactions that occurred (from gillespie)
    xs = xs[:, 1:(end-1)] .* reshape(dts, 1, length(dts)) #scale xs by time
    coords = 1000:1000:length(dts) #timepoints to plot; plot every 1000 iterations
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
    figure(figsize = (5,3.5)) #plot fluxes
    cols = ["r", "b", "g"]
    for i in 1:3 #x1, x2 and x3 on the same graph
        plot(coords, normfluxes[i,:], cols[i]*"--")
    end
    plot([0; coords[end]], [0; 0], "k--")
    legend(["x1", "x2", "x3"])
    ylabel(L"\dfrac{<R_x^+> - <R_x^- >}{<x>}")
    xlabel("iteration")
    ylim(ylimvals[1], ylimvals[2])
    savefig(fname, bbox_inches="tight", dpi=120)
    close()
end

function get_noise(dts, xs; tau1=1, tau2=2, tau3=4, alpha=10, lambda=1, Print=true)
    #return intrinsic and extrinsic noise for the 3 components
    #return sigma_x^2 / <x>^2 for the three components
    #return the full covariance matrix
    #return the expected analytical covariance matrix
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
        intrinsic[i] = 1/meanx #intrinsic noise (poisson)
        if i == 1
            expec = 1/meanx #only poisson
            extrinsic[i] = 0 #no extrinsic noise
        elseif i == 2
            #calculate expected value from part B
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
    anal_eta = zeros(3,3) #calculate expected analytical values
    anal_eta[1,1] = 1/(alpha*tau1)
    anal_eta[2,2] = 1/(alpha*lambda*tau1*tau2) + 1/(alpha*(tau1+tau2))
    anal_eta[3,3] = 1/(alpha*lambda*tau1*tau3) + 1/(alpha*(tau1+tau3))
    anal_eta[1,2] = 1/(alpha*(tau1+tau2)); anal_eta[2,1] = anal_eta[1,2]
    anal_eta[1,3] = 1/(alpha*(tau1+tau3)); anal_eta[3,1] = anal_eta[1,3]
    anal_eta[2,3] = tau3/(alpha*(tau1+tau3)*(tau2+tau3))+tau2/(alpha*(tau1+tau2)
                    *(tau3+tau2));   anal_eta[3,2] = anal_eta[2,3]

    return intrinsic, extrinsic, noises, eta, anal_eta #return var/mean^2
end

function plot_figs(ts, xs, xweights; fname="figures/test")
    #plot part of trajectory and histogram of x values
    figure() #plot part of the trajectory
    t1, t2 = 1000, 1500
    plot(ts[t1:t2], xs[1,t1:t2], "r--") #x1
    plot(ts[t1:t2], xs[2,t1:t2], "b--") #x2
    plot(ts[t1:t2], xs[3,t1:t2], "g--") #x3
    savefig(fname*"_trajec.png")
    close()

    figure() #plot a histogram
    cols = ["r", "b", "g"]
    for i = 1:3
        plot(0:(size(xweights)[1]-1), xweights[:,i], cols[i]*"-")
    end
    legend(["x1", "x2", "x3"])
    xlim(0, 40)
    savefig(fname*"_hist.png")
    close()
end

function plot_param_results(fname, param, values, ints, exts, errs,
                            var_errors, CV32s, rho23s, all_noises)
    #plot graphs for questions 2A, 2B, 2C and 2E

    #plot intrinsic and extrinsic noise against parameter values
    figure(figsize = (5,3.5))
    cols = ["r", "b", "g"]
    for i in 2:3 #for x2 and x3
        plot(values, ints[i,:], cols[i]*"-")
        plot(values, exts[i,:], cols[i]*"--")
    end
    legend(["x2 intrinsic", "x2 extrinsic",
            "x3 intrinsic", "x3 extrinsic"])
    xscale("log"); yscale("log")
    xlabel(param); ylabel("normalized variance")
    savefig(fname*"_int_ext.png", bbox_inches="tight", dpi=120)
    close()

    #plot ratio of intrinsic to extrinsic noise
    figure(figsize = (5,3.5))
    for i in 2:3 #for x2 and x3
        plot(values, ints[i,:]./exts[i,:], cols[i]*"-")
    end
    legend(["x2 intrinsic/extrinsic",
            "x3 intrinsic/extrinsic"])
    xscale("log"); yscale("log")
    xlabel(param); ylabel("variance")
    savefig(fname*"_ratio.png", bbox_inches="tight", dpi=120)
    close()

    #plot scatterplot of intrinsic vs extrinsic noise
    for i=2:3 #separate plots for x2 and x3
        figure(figsize = (5,3.5))
        scatter(x=ints[i,:], y=exts[i,:], c=["" "b" "g"][i], s=2)
        xlabel("Intrinsic")
        ylabel("Extrinsic")
        title("X"*string(i)*" noise, varying lambda")
        xscale("log")
        yscale("log")
        xlim(0.7*minimum(ints[i,:]), 1.5*maximum(ints[i,:]))
        if param == "lambda" ylim(10^-2,3*10^-2)
        else ylim(0.7*minimum(exts[i,:]), 1.5*maximum(exts[i,:])) end
        savefig(fname*"_scatter_x"*string(i)*".png",bbox_inches="tight", dpi=120)
    end

    #plot observed vs. expected noise
    for i=1:3 #both x1 x2 and x3
        figure(figsize = (5,3.5))
        scatter(x=ints[i,:].+exts[i,:],y=all_noises[i,:],c=["r" "b" "g"][i],s=2)
        xlabel(L"$\eta_{kk}^{exp}$")
        ylabel(L"$\eta_{kk}^{obs}$")
        xscale("log")
        yscale("log")
        xlim(0.7*minimum(ints[i,:].+exts[i,:]), 1.5*maximum(ints[i,:].+exts[i,:]))
        ylim(0.7*minimum(all_noises[i,:]), 1.5*maximum(all_noises[i,:]))
        savefig(fname*"_Q2C_x"*string(i)*".png",bbox_inches="tight", dpi=120)
    end

    #plot histograms of flux errors
    for i in 1:3 #for x1, x2 and x3
        figure(figsize = (5, 3.5))
        PyPlot.plt[:hist](errs[i,:], 40, color=cols[i])
        xlabel("relative flux error")
        ylabel("frequency")
        maxx = maximum(abs.(errs[i,:]))*1.05
        xlim(-maxx, maxx)
        savefig(fname*"_hist"*string(i)*".png", bbox_inches="tight", dpi=120)
        close()
    end

    #plot histograms of deviations of covariances from analytical expressions
    inds = [[1 1],[2 2],[3 3],[1 2],[1 3],[2 3]] #matrix indices
    for i in 1:6 #six histograms
        figure(figsize = (5, 3.5))
        PyPlot.plt[:hist](var_errors[i,:], 20)
        lab = "eta"*string(inds[i][1])*string(inds[i][2])
        xlabel(lab)
        ylabel("frequency")
        maxx = maximum(abs.(var_errors[i,:]))*1.05
        xlim(-maxx, maxx)
        savefig(fname*"_hist_eta"*string(i)*".png", bbox_inches="tight", dpi=120)
        close()
    end

    #plot rho23 vs CV3/CV2
    figure(figsize = (5,3.5))
    plot(CV32s, rho23s, "kx")
    xlabel(L"CV3/CV2"); ylabel(L"\rho_{23}")
    savefig(fname*"_rho23_CV32.png", bbox_inches="tight", dpi=120)
    close()
end

function vary_params(param, values; fname="figures/test_noise",
                    tau1=3, tau2=2, tau3=4, lambda=1, alpha=10, niter=500000)
    #calculate & plot intrinsic and extrinsic noise for a range of parameters
    #param is name of parameter to vary, values is set of values
    N = length(values) #number of parameter values
    exts = zeros(3, N); ints = zeros(3, N) #store extrinsic and intrinsic noise
    errs = zeros(3,N) #store relative flux errors
    var_errors = zeros(6,N) #store n11, n22, n33, n12, n13, n23
    inds = [[1 1],[2 2],[3 3],[1 2],[1 3],[2 3]] #corresponding matrix indices
    rho23s = zeros(N) #pearson correlations for x2 and x3
    CV32s = zeros(N) #CV3/CV2
    all_noises = zeros(3,N) #variances
    t = time()
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
        intrinsic, extrinsic, noises, eta, anal_eta =
                    get_noise(dts, xs; tau1=tau1, tau2=tau2, tau3=tau3,
                            alpha=alpha, lambda=lambda, Print=false)

        #store intrinsic and extrinsic noise
        for j in 1:3
            ints[j, i] = intrinsic[j]
            exts[j, i] = extrinsic[j]
            #also calculate flux deviation from 0
            errs[j, i] = (sum(rxns[2*j-1,:]) - sum(rxns[2*j,:])
                        ) / sum(xs[j, 1:(end-1)] .* dts)
        end

        #calculate relative covariance errors
        for j in 1:6
            var_errors[j, i] = (eta[inds[j][1], inds[j][2]] -
                            anal_eta[inds[j][1], inds[j][2]]) /
                            anal_eta[inds[j][1], inds[j][2]]
        end

        #calculate rho23 and CV3/CV2
        rho23s[i] = cor(xs[2,:], xs[3,:])
        CV32s[i] = sqrt(noises[3]/noises[2]) #CV2 is sqrt(eta_22)
        all_noises[:,i] = noises #also store CV1^2, CV2^2 and CV3^2

        #print progress
        println(param, " ", i, " ", val, "  time:", time()-t, "  err: ", errs[:,i])
        t = time()
    end

    names = ["values" "ints" "exts" "errs" "var_errors" "CV32s" "rho23s" "noises"]
    mats = [values, ints, exts, errs, var_errors, CV32s, rho23s, all_noises]
    for i in 1:length(names) #save data for future replotting
        writedlm("dlms/"*fname*"_"*names[i]*".dlm", mats[i])
    end
    #plot all of our results
    plot_param_results("figures/"*fname, param, values, ints, exts, errs,
                        var_errors, CV32s, rho23s, all_noises)
end

function plot_exp_rho_CVs(lambdas = 10 .^ (-2:(4/99):2))
    #plot analytical values of rho23 and CV3/CV2
    N = length(lambdas)
    rhos = zeros(N)
    CVs = zeros(N)
    for i = 1:N #for each value of lambda
        l = lambdas[i]
        rhos[i] = 17/105*sqrt(2520*l^2/(72*l^2+102*l+42))
        CVs[i] = sqrt(120*l/(168*l+140) + 70/(168*l+140))
    end
    figure(figsize = (5,3.5)) #plot result
    plot(lambdas, rhos, "b--")
    plot(lambdas, CVs, "g--")
    xscale("log")
    legend([L"\rho_{23}", L"CV3/CV2"])
    xlabel(L"\lambda")
    savefig("figures/exp_rho_CV.png", dpi=120, bbox_inches="tight")
    close()
end

#run actual simulations and plot figures


println("\nnew simulation")

#step sizes
step_size = [ [1; 0; 0], [-1; 0; 0], #x1 +- 1
            [0; 1; 0], [0; -1; 0], #x2 +- 1
            [0; 0; 1], [0; 0; -1]] #x3 +- 1
tau2=2 #fixed
tau3=4 #fixed
#variable parameters
tau1=3
lambda=1
alpha=10
niter = 5*10^7 #number of iterations

indrun = false
if indrun #run simulations and plot timecourses
    fnames = ["01" "100"]
    vals = [0.01 100] #lambda values
    ylimvals = [[-0.02 0.02], [-0.3 0.3]]
    for i = 1:2 #for each value of lambda
        #start with equilibrium conditions for improved convergence
        x0 = [Int(round(alpha*tau1));
            Int(round(alpha*tau1*tau2*vals[i]));
            Int(round(alpha*tau1*tau3*vals[i]))]
        #run simulation
        dts, ts, xs, xweights, rxns = gillespie(x0, step_size, niter=500000,
                                    alpha=alpha, lambda=vals[i],
                                    tau1=tau1, tau2=tau2, tau3=tau3)
        #calculate noise
        intrinsic, extrinsic, noises, eta, anal_eta =
                get_noise(dts, xs, tau1=tau1, tau2=tau2, tau3=tau3,
                        alpha=alpha, lambda=vals[i])
        #plot timecourse of concentrations and flux balance
        plot_figs(ts, xs, xweights, fname="figures/test_lambda"*fnames[i]*".png")
        plot_flux(rxns, xs, dts, ts, fname="figures/flux_lambda"*fnames[i]*".png",
                    ylimvals = ylimvals[i])
    end
end

#run simulations for different value of lambda
#vary_params("lambda", 10 .^ (-2:(4/99):2),
#            fname="figures/iter"*string(niter)*"_lambda_noise",
#            niter=niter)

#find noise vs. lambda
vary_params("lambda", 10 .^ (-2:(4/99):2),
            fname="iter"*string(niter)*"_lambda2",
            niter=niter )

#run simulations for different value of alpha
vary_params("alpha", 10 .^ (-2:(3/99):1),
            fname="iter"*string(niter)*"_alpha",
            niter=niter)

#run simulations for different value of tau1
vary_params("tau1", 10 .^ (-2:(5/99):3),
            fname="iter"*string(niter)*"_tau1",
            niter=niter)

plot_exp_rho_CVs() #plot analytical values of rho23 and CV3/CV2
