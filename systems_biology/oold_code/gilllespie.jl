#Gillespie??

using PyCall, PyPlot, LinearAlgebra, Distributions,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats

function rates_LV(x)
    #x = x, y
    alpha = 0.1; beta = 0.1
    delta = 0.1; gamma = 0.1


    return [alpha*x[1]
            beta*x[1]*x[2]
            delta*x[1]*x[2]
            gamma*x[2]]

end

function gillespie(x0, deltas; tmax = 1000)
    #x0, deltas are colvectors
    t = 0
    x = x0 #colvector
    N = length(deltas)
    xs = zeros(length(x), 1); xs[:,1] = x0
    ts = zeros(1, 1); ts[1,1] = t
    while (t < tmax)
        rates = rates_LV(x) #calculate new rates
        rtot = sum(rates) #total rate
        dt = rand(Exponential(rtot)) #time till next event
        i = sample(1:N, Weights(rates ./ rtot)) #draw reaction to occur
        t += dt #update time
        x += deltas[i] #update state
        ts = [ts t] #store new time
        xs = [xs x] #store new state
    end
    return ts, xs
end

x0 = [10; 8]
deltas = [ [1; 0], [-1; 0], [0; 1], [0; -1] ]
ts, xs = gillespie(x0, deltas, tmax=10000)

figure()
plot(ts[1,:], xs[1,:], "r-")
plot(ts[1,:], xs[2,:], "b-")
savefig("test.png")
