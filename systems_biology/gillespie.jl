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
    xcounts = zeros(1000,Nx)
    for k in 1:Nx xcounts[x[k]+1, k] += 1 end
    #ts = zeros(1, 1); ts[1,1] = t
    #while (t < tmax)
    for n in 1:niter
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
    return dts, ts, xs, xcounts
end

function get_noise(dts, xs; tau1=1, tau2=2, tau3=4)
    #return sigma_x / <x>
    N = size(xs)[1]
    x2s = xs.^2
    ttot = sum(dts)
    noises = zeros(N)
    means = zeros(N)
    for i in 1:3
        meanx = sum(xs[i, 1:(end-1)] .* reshape(dts, length(dts), 1))/ttot
        meanx2 = sum(x2s[i, 1:(end-1)] .* reshape(dts, length(dts), 1))/ttot
        println(meanx2," ", meanx)
        varx = meanx2 - meanx^2
        normstd = sqrt(varx)/meanx
        #normstd = std(xs[i,:])/meanx
        noises[i] = normstd
        means[i] = meanx
        if i == 1
            expec = 1/meanx
        elseif i == 2
            expec = 1/meanx+tau1/((tau1+tau2)*means[1])
        elseif i == 3
            expec = 1/meanx+tau1/((tau1+tau3)*means[1])
        else
            println("something is wrong...")
        end
        println(i, "  normstd: ", normstd, "  exp: ", sqrt(expec))
    end


    return noises
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


dts, ts, xs, xcounts = gillespie(x0, step_size, niter=100000, alpha=1,
                            tau1=tau1, tau2=tau2, tau3=tau3)

print(get_noise(dts, xs, tau1=tau1, tau2=tau2, tau3=tau3))

figure()
t1, t2 = 1000, 1500
plot(ts[t1:t2], xs[1,t1:t2], "r--")
plot(ts[t1:t2], xs[2,t1:t2], "g--")
plot(ts[t1:t2], xs[3,t1:t2], "b--")
savefig("test_trajec.png")
close()

xcounts /= sum(xcounts[:,1])
xvals = 0:20
figure()
cols = ["r", "g", "b"]
lambdas = [10 5 3]
betas = [1 2 3]
for i = 1:3
    lambda = lambdas[i]; beta = betas[i]
    plot(0:(size(xcounts)[1]-1), xcounts[:,i], cols[i]*"-")
    #ps = (lambda/beta).^xvals .* exp.(-lambda/beta) ./ (factorial.(xvals))
    #plot(xvals, ps, cols[i]*"--")
end
legend(["x1", "x2", "x3"])
xlim(0, 20)
savefig("test_hist.png")
close()
