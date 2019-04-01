#Code for simulatign and analyzing occular dominance

using PyCall, PyPlot, LinearAlgebra, Distributions, FFTW,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats

function plot_occ(wm; xlab="", ylab="", filename="default", Title="default")
        #plots the degree of occular dominance as a 1d heatmap
        Wm = zeros(1,512)
        Wm[1,:] = wm #1x512 array
        figure(figsize = (10,1.5))
        imshow(Wm, cmap=ColorMap("gray_r"), vmin = -1, vmax = +1, aspect="auto")
        print("plotted")
        xticks([], [])
        yticks([], [])
        savefig(filename, bbox_inches = "tight")
        close()
end

function get_int(d1, d2; n = 512, sigma = 0.066, L = 10)
    #return K(|d1-d2|) as specified in the question
    Kij = exp( -( d1 - d2 )^2/(2*sigma^2) ) -
            1/9 * exp( -( d1 - d2 )^2/(18*sigma^2) )
    return Kij
end

function get_trans(k; n = 512, sigma = 0.066, L = 10)
    #for a given value of k, returns the fourier transform of K for this k
    Kij = sqrt(2*pi*sigma^2)*exp(-k^2*sigma^2/2) -
        1/9 * sqrt(18*pi*sigma^2) * exp(-9*k^2*sigma^2/2)
    return Kij
end

function plot_Ktildes()
    #plots the fourier transform of K as a function of k for a range
    #of mu values
    mus = 0:1:96
    ks = [2*pi/10*mu for mu in mus]
    Ks = [get_trans(k) for k in ks]
    figure(figsize = (5,3))
    plot(ks, Ks)
    xlabel("k (1/mm)")
    ylabel(L"\tilde K")
    xlim(minimum(ks), maximum(ks))
    ylim(0, maximum(Ks*1.05))
    savefig("figures/plot_Ktilde.png", bbox_inches = "tight")
    close()
    return ks, Ks
end

function plot_Ks()
    #plots K(0, d) for a range of intercortical distances
    ds = -0.6:0.01:0.6
    Ks = [get_int(d, 0) for d in ds]
    figure(figsize = (5,3))
    plot(ds, Ks)
    xlabel("cortical distance (mm)")
    ylabel("K")
    xlim(minimum(ds), maximum(ds))
    ylim(minimum(Ks)-0.05, maximum(Ks)+0.05)
    savefig("figures/plot_K.png", bbox_inches = "tight")
    close()
    return ds, Ks
end

function plot_eigvec(;mu=13)
    #plots the eigenvector corresponding to a given eigenvalue. 13 is
    #the principal eigenvalue
    ds = -0.6:0.01:0.6
    vals = [cos(2*pi*mu* d*512/10 /512) for d in ds]
    figure(figsize = (5,3))
    plot(ds, vals)
    xlabel("cortical distance (mm)")
    ylabel(L"e")
    xlim(minimum(ds), maximum(ds))
    ylim(-1.05, 1.05)
    savefig("figures/plot_eigvec.png", bbox_inches = "tight")
    close()
    return ds, vals
end

function update(W; epsilon = 0.01, K=K, Q=Q)
    #updates the weights a single iteration with learning rate epsilon
    #K: 512x512, W: 512x2, Q: 2x2
    W += epsilon * K * W * Q
    W .+= (1 .- sum(W, dims=2))*0.5 #explicit normalization step
    W[ W.<0 ] .= 0 #threshold at zero
    W[ W.>1 ] .= 1 #threshold at 1
    return W
end

function run_sim(;n = 512, sigma = 0.066, L = 10, N=1000, temp = true,
                fname="occsimtemp", Plot=true, q_D=0.7)
    #runs simulation of cortical columns for N timesteps.
    K = zeros(n, n)
    for i in 0:(n-1)
        for j in 0:(n-1)
            #construct cortical interaction matrix
            K[i+1, j+1] = get_int( L/(n-1)*i, L/(n-1)*j, sigma=sigma )
        end
    end
    d = Normal(0, 0.01) #gaussian noise
    wLs = 0.5 .+ rand(d, 512) #initialize left weights
    wRs = 1.0 .- wLs #sum of weights is always 1.
    W = [wLs wRs]
    Q = ones(2,2); Q[1,2] = q_D; Q[2,1] = q_D #construct correlation matrix
    stds = zeros(N) #use std as proxy for convergence
    for i = 1:N
        W = update(W, K=K, Q=Q) #update weights
        Wm = W[:,2]-W[:,1] #calculate w_(-)
        stds[i] = std(Wm)
        if temp
            if i%50 == 0 #plot w_(-) every 50 iterations
                plot_occ(Wm, filename="figures/"*fname*"_int/i"*string(i)*".png")
            end
        end
    end

    Wm = W[:,2]-W[:,1]
    Plot && plot_occ(Wm, filename="figures/"*fname*"_heat.png") #plot result

    figure(figsize=(5,3)) #plot standard deviations
    plot(1:N, stds)
    xlabel("iteration")
    ylabel(L"std({\bf w_-})")
    ylim(0,1)
    savefig("figures/"*fname*"_stds.png", bbox_inches="tight")
    close()

    return W, Wm
end

function plot_mags(mags)
    #plot the magnitudes of the discrete fourier transform of w_(-)
    mus = 0:1:96
    ks = [2*pi/10*mu for mu in mus]
    Ks = [mags[mu+1] for mu in mus]
    figure(figsize = (5,3))
    plot(ks, Ks)
    xlabel("k (1/mm)")
    ylabel("|DFT|")
    xlim(minimum(ks), maximum(ks))
    ylim(0, maximum(Ks*1.05))
    savefig("figures/plot_DFT.png", bbox_inches = "tight")
    close()

end

function repeat_sims(N = 5000)
    #run N simulations, calculate the mean DFT magnitudes and plot this

    dfts = zeros(512,N)
    for i = 1:N
        println("new i: ",i)
        W, Wm = run_sim(Plot=false, temp=false)
        dft = fft(Wm) #perform a fast fourier transform
        mags = norm.(dft) #get magnitude of components
        dfts[:,i] = mags
    end

    mags_mean = mean(dfts, dims=2) #get mean
    plot_mags(mags_mean) #plot result
    writedlm("mags_mean.dlm", mags_mean) #write result to file
    #print optimum value of mu
    println((0:96)[mags_mean[1:97] .== maximum(mags_mean[1:97])], " : ",  maximum(mags_mean[1:97]) )
    return mags_mean
end

function get_max_k(sig)
    #for a given value of sigma, finds the value of k that maximizes K_tilde
    ks = 0:0.01:100
    Ks = [get_trans(k, sigma=sig) for k in ks] #get fourier transform at k
    kmax = ks[ Ks .== maximum(Ks) ][1] #find max
    return kmax
end

function scan_widths(sigs = 0.01:0.001:0.200)
    #for a range of sigma values, finds the k value that optimizes
    #K_tilde and plots this result
    N = length(sigs)
    ks = zeros(N)
    for i = 1:N
        k = get_max_k(sigs[i]) # find max k value
        ks[i] = k
    end
    figure(figsize = (5,3)) #plot result
    plot(sigs, ks)
    xlim(minimum(sigs), maximum(sigs))
    ylim(0, maximum(ks))
    plot([0, 1], [8.17, 8.17], "b:")
    ylabel(L"k_{max}")
    xlabel(L"\sigma")
    savefig("figures/test_maxk.png", bbox_inches="tight")
    close()
    return sigs, ks
end

W, Wm = run_sim(fname="occsim_066", sigma=0.066, N=3000) #default parameters
sigs, ks = scan_widths() #find how stripe width varies with k
W, Wm = run_sim(fname="occsim_012", sigma=0.012, N=3000) #narrow stripes
W, Wm = run_sim(fname="occsim_174", sigma=0.174, N=3000) #wide stripes

mags = repeat_sims() #run repeated simulation to get |DFT|

#Plot K, it's fourier transform, and it's principal eigenvector
ds, Ks = plot_Ks()
ks, Ktildes = plot_Ktildes()
ds, vals = plot_eigvec()
