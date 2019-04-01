##code for the cart pole balancing problem
using PyCall, PyPlot, LinearAlgebra, Distributions,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats

#Physical parameters
#g = -9.8 #m/s2, acceleration due to gravity
g = 9.8
mc = 1.0 #kg, mass of cart
mp = 0.1 #kg, mass of pole,
l = 0.5 #m, half-pole length,
muc = 0.0005#, coefficient of friction of cart on track,
mup = 0.000002#, coefficient of friction of pole on cart,
Ft = 10.0 # p/m 10 newtons, force applied to cart's center of mass at time t.

#model parameters for Barto et al.
alpha = 1000
beta = 0.5
delta = 0.9
gamma = 0.95
lambda = 0.80
sigma = 0.01
dnorm = Normal(0, sigma^2)
dt = 0.02

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


function plot_poles(N, xs, thetas, performance, dt, fname, vs; all=true)
    #N is number of timesteps, xs and thetas are given for a simulation
    #performance is simulation length by trial number

    if all #make all plots
    #plot cart position vs time
    ts = (1:N)*dt
    figure(figsize = (5,3))
    plot(ts, xs)
    plot(ts, repeat([-2.4], N), "k--")
    plot(ts, repeat([2.4], N), "k--")
    xlabel("time (s)")
    ylabel("x")
    savefig("figures/"*fname*"_xs.png", bbox_inches = "tight")
    close()

    #plot theta vs time
    figure(figsize = (5,3))
    plot(ts, thetas)
    plot(ts, repeat([-12], N), "k--")
    plot(ts, repeat([12], N), "k--")
    xlabel("time (s)")
    ylabel("theta")
    savefig("figures/"*fname*"_thetas.png", bbox_inches = "tight")
    close()

    #plot both x and theta
    figure(figsize = (5,3))
    plot(ts, xs*(1/2.4), "b:")
    plot(ts, thetas*(1/12), "r:")
    plot(ts, repeat([-1], N), "k--")
    plot(ts, repeat([1], N), "k--")
    legend(["cart position", "pole angle"])
    xlabel("time (s)")
    ylabel("theta")
    savefig("figures/"*fname*"_xs_thetas.png", bbox_inches = "tight")
    close()

    if length(vs) > 5
        #plot heatmap of predicted reward
        newvs = zeros(3, 6)
        for (i,thetadot) = enumerate([-100;0;100])
            for (j,theta) = enumerate([-9; -3.5; -0.5; 0.5; 3.5; 9])
                temps = zeros()
                for x = [-1.6;0;1.6] #average over x and xdot
                    for xdot = [-1;0;1]
                        temps = [temps; vs[get_state([x;theta;xdot;thetadot])]]
                    end
                end
                newvs[i,j] = mean(temps)
            end
        end
        heatmap(newvs, [-9; -3.5; -0.5; 0.5; 3.5; 9], [-100;0;100];
                xlab=L"$\theta$", ylab=L"$\dot \theta$",
                filename="figures/"*fname*"_vs.png", Title="Predicted Reward")
    end

    end

    #plot performance vs trial number
    ntrial = length(performance)
    figure(figsize = (5,3))
    plot(1:ntrial, performance)
    xlabel("trial")
    ylabel("time (s)")
    savefig("figures/"*fname*"_performance.png", bbox_inches = "tight")
    close()
end

function smooth_performance(performance, dt)
    #compute rolling average of trial length as in Barto et al.
    N = length(performance)
    newperf = zeros(N)
    for i = 1:4
        newperf[i] = mean(performance[1:i])
    end
    for i = 5:N #average over 5 timepoints
        newperf[i] = mean(performance[i-4:i])
    end
    return newperf*dt #real time
end

function get_state(state)
    #return unique state index given vector of [x,theta,xdot,thetadot]
    x, theta, xdot, thetadot = state

    k1 = sum(x .> [-0.8, 0.8])
    k2 = sum(theta .> [-6,-1,0,1,6])
    k3 = sum(xdot .> [-0.5, 0.5])
    k4 = sum(thetadot .> [-50, 50])
    k = 54*k1 + 9*k2 + 3*k3 + k4 + 1
    return k #index of non-zero element
end

function get_output(state, ws; d=dnorm)
    #get control action. +1 is right, -1 is left
    k = get_state(state)
    y = sign(ws[k] + rand(d)) #get output
    return y, k
end

function update_es(es, k, y, delta)
    #Update ASE eligibility trace. k is index of current state
    es *= delta #decay
    es[k] += (1-delta)*y #eligibility trace increased for current state
    return es
end

function update_xbars(xbars, k, lambda)
    #Update ACE eligibility trace. k is index of current state.
    xbars *= lambda #decay
    xbars[k] += (1-lambda) #eligibility trace increased for current state
    return xbars
end

function update_ws(ws, r, es, alpha)
    #update ASE weights
    ws += alpha*r*es
    return ws
end

function update_vs(vs, r, xbars, beta)
    #update ACE weights
    vs += beta*( r )*xbars
    return vs
end

function update_rhat(r, pt, ptm, gamma)
    #pt = p(t). ptm=p(t-1)
    rhatt = r + gamma*pt-ptm
    return rhatt
end

function update_cart(state, yt; dt = 0.02, l=l, mc=mc, mp=mp, muc=muc, mup=mup, g=g)
    #update the physical system according to equations in Barto et al.'s appendix
    x, theta, xdot, thetadot = state
    theta, thetadot = 2*pi/360*theta, 2*pi/360*thetadot #convert to radians
    Ft = 10*yt #force to apply

    dthetadot = ( g*sin(theta) +
                cos(theta)*( ( -Ft-mp*l*thetadot^2*sin(theta) + muc*sign(xdot) )/( mc + mp ) ) -
                mup*thetadot/(mp*l) ) / (
                l*( 4/3 - ( mp* (cos(theta)^2) )/( mc+mp ) ))
    dxdot = (  Ft + mp*l*(thetadot^2*sin(theta) - dthetadot*cos(theta)) -
            muc*sign(xdot) ) / ( mc+mp )

    #update state
    thetadot += dthetadot*dt
    xdot += dxdot*dt
    theta += thetadot*dt
    x += xdot*dt
    return [x; 360/(2*pi)*theta; xdot; 360/(2*pi)*thetadot]
end

function state_failed(state)
    #find out if a given state is outside the system boundaries
    if abs(state[1]) > 2.4 || abs(state[2]) > 12
        return true #system is in a failure state
    end
    return false #state is allowed

end

function update_2layer(state, ws, vs, fs, D, cs, A, ptm, oldfailed; d=dnorm,
                        alpha=alpha, lambda=lambda, gamma=gamma, delta=delta, beta=beta,
                        l=l, mc=mc, mp=mp, muc=muc, mup=mup, layer1=false)
    #perform a single iteration of the 2-layer system from Anderson 1987
    #if 1layer = true, use only a single layer of this system

    #get input activities. x5 is constant bias term
    xs = [(state[1]+2.4)/4.8;(state[2]+12)/24;(state[3]+1.5)/3;(state[4]+115)/230;0.5]
    zs = 1 ./ (1 .+ exp.(-D*xs)) #ASE hidden layer

    P = ws' * xs + fs' * zs
    P = 1 / (1 + exp(-P)) #probability of rightwards impulse
    if rand() < P
        yt = 1
    else
        yt = -1
    end

    qstt = 1 ./ (1 .+ exp.(-A*xs)) #ACE hidden layer
    ptt = vs' * xs + cs' * qstt #predicted reward

    state = update_cart(state, yt, l=l, mc=mc, mp=mp, muc=muc, mup=mup) #state t+1
    #get new input values
    newxs = [(state[1]+2.4)/4.8;(state[2]+12)/24;(state[3]+1.5)/3;(state[4]+115)/230;0.5]

    qstp = 1 ./ (1 .+ exp.(-A*newxs)) #q[t,t+1]
    ptp = vs' * newxs + cs' * qstp #p[t,t+1]

    #calculated ACE reward signal
    if oldfailed #we failed on our previous iteration
        rhatt = 0
        failed = false
    elseif state_failed(state) #current state is a failed state
        rhatt = -1 - ptt
        failed = true
    else
        rhatt = gamma*ptp - ptt
        failed = false
    end
    qs = qstt
    pt = ptt

    #update ASE weights
    ws += alpha*rhatt*(max(yt,0) - P)*xs
    D += ( ((0.2*alpha) * rhatt * zs) .* (1 .- zs) .* sign.(fs) * (max(yt,0) - P) ) * xs'
    fs += alpha * rhatt * (max(yt,0) - P)*zs

    #update ACE weights
    vs += beta*rhatt*xs
    A += ( ((0.25*beta)*rhatt*qs) .* (1 .- qs) .* sign.(cs) ) * xs'
    cs += beta * rhatt * qs

    if failed #initialize new state
        state = [rand()*4.8-2.4;rand()*24-12;0;0]
    end
    if layer1 #set second layer weights to 0
        D = zeros(5,5); fs = zeros(5); A = zeros(5,5); cs = zeros(5)
    end

    return state, ws, vs, fs, D, cs, A, pt, yt, P, zs, qs, xs, failed
end

function update_system(state, ws, vs, es, xbars, ptm; d=dnorm,
                        alpha=alpha, lambda=lambda, gamma=gamma, delta=delta, beta=beta,
                        l=l, mc=mc, mp=mp, muc=muc, mup=mup)
    #update the system from Barto et al.

    yt, k = get_output(state, ws) #get impulse at time t
    pt = vs[k] #get probability of future reward at time t
    if state == [0;0;0;0] rhatt = 0 else
    rhatt = gamma*pt - ptm #calculate expected reward at time t
    end

    ws = update_ws(ws, rhatt, es, alpha) #update ASE weights t+1
    vs = update_vs(vs, rhatt, xbars, beta) #update ACE weights t+1
    es = update_es(es, k, yt, delta) #update ASE eligibility traces t+1
    xbars = update_xbars(xbars, k, lambda) #update ACE eligibility traces t+1

    state = update_cart(state, yt, l=l, mc=mc, mp=mp, muc=muc, mup=mup) #update states to t+1

    return state, ws, vs, es, xbars, pt, k
end

function run_sim(;fname="test", N = 500000, ntrials = 100, mc=mc, mp=mp, l=l, muc=muc,
                mup=mup, layer2 = false, alpha=alpha, beta=beta, gamma=gamma,
                lambda=lambda, delta=delta, d=dnorm, Plot=true, Print=true,
                all = false, layer1 = false)
    #run a simulation for N timesteps of ntrials trials, whichever occurs first
    #layer2=true for an Anderson simulation, layer1=true to do this with 1 layer

    n = 162 #number of Barto states
    pt = 0 #no initial predicted reward

    if layer2 #Anderson
        A, D = zeros(5,5), zeros(5,5)
        ws, vs, fs, cs = zeros(5), zeros(5), zeros(5), zeros(5)
    else
        ws, vs, es, xbars = zeros(n), zeros(n), zeros(n), zeros(n)
    end
    state = [0;0;0;0] #start upright
    ntrial, lasttrial, performance = 0, 0, zeros(ntrials) #keep track of performance
    thetas = zeros(N) #store angles
    xs = zeros(N) #store positions
    Print && println("\nnew simulation")
    failed = true #start from a newly initialized state; history irrelevant
    for i = 1:N
        if !layer2 #Barto update step
        state, ws, vs, es, xbars, pt, k = update_system(state, ws, vs, es, xbars, pt,
                                                    l=l, mc=mc, mp=mp, muc=muc, mup=mup,
                                                    alpha=alpha, beta=beta, gamma=gamma,
                                                    lambda=lambda, delta=delta)
        else #Anderson update step
        state, ws, vs, fs, D, cs, A, pt, yt, P, zs, qs, xvals, failed = update_2layer(state, ws, vs, fs, D, cs, A, pt, failed,
                                                        l=l, mc=mc, mp=mp, muc=muc, mup=mup,
                                                        alpha=alpha, beta=beta, gamma=gamma,
                                                        lambda=lambda, delta=delta, layer1 = layer1)
        end
        xs[i] = state[1] #store data
        thetas[i] = state[2]
        if failed #if failed Anderson state, update trial number
            ntrial += 1
            Print && println("system failed i=",i,"  ntrial=",ntrial)
            performance[ntrial] = i-lasttrial
            lasttrial = i
            if ntrial == ntrials break end
        end
        if (abs(state[1]) >= 2.4 || abs(state[2]) >= 12) && !layer2 #failed Barto
            ntrial += 1
            #print what happened
            if state[1] <= -2.4
                Print && println("cart hit left, i=",i,"  ntrial=",ntrial)
            elseif state[1] >= 2.4
                Print && println("cart hit right, i=",i,"  ntrial=",ntrial)
            elseif state[2] <= -12
                Print && println("pole fell left, i=",i,"  ntrial=",ntrial)
            else
                Print && println("pole fell right, i=",i,"  ntrial=",ntrial)
            end

            #we're now not in any of the 162 states
            rhatt = -1-pt
            pt = 0
            #learn from error signal
            ws = update_ws(ws, rhatt, es, alpha) #update ASE weights t+1
            vs = update_vs(vs, rhatt, xbars, beta) #update ACE weights t+1
            es = zeros(n) #update_es(es, k, yt, delta) #update ASE eligibility traces t+1
            xbars = update_xbars(xbars, k, lambda)

            state = [0;0;0;0] #-1 #reset system
            performance[ntrial] = i-lasttrial #store performance of last trial
            lasttrial = i
            if ntrial == ntrials break end
        end
    end

    if ntrial < ntrials #if we reached iteration threshold, extrapolate performance
        performance[ntrial+1:end] .= max(performance[ntrial], N-lasttrial)
    end
    performance = smooth_performance(performance, dt) #smooth performance as in Barto
    Plot && plot_poles(N, xs, thetas, performance, dt, fname, vs; all=all)

    return performance
end

function test_performance(;n=10, ntrials=50, mc=mc, mp=mp, l=l, muc=muc, mup=mup,
                fname="test", Plot=true, alpha=alpha, beta=beta, gamma=gamma,
                lambda=lambda, delta=delta, d=dnorm)
    #run n simulations with maximum ntrial trials and average performance
    perfs = zeros(ntrials, n) #performances
    learns = zeros(n) #trials needed to learn
    for i = 1:n
        perf = run_sim(;ntrials = ntrials, mc=mc, mp=mp, l=l, muc=muc, mup=mup,Plot=false,
                alpha=alpha, beta=beta, gamma=gamma, lambda=lambda, delta=delta, Print=true)
                #run simulation
        perfs[:,i] = perf
        learn = sum(perf .< 500000/10*0.02) #define learning as threshold performance
        learns[i] = learn
    end

    avgperf = mean(perfs, dims=2)
    if Plot #plot average performance
        figure(figsize = (5,3))
        plot(1:ntrials, avgperf)
        xlabel("trial")
        ylabel("time (s)")
        savefig("figures/"*fname*"_meanperformance.png", bbox_inches = "tight")
        close()
    end

    learn = mean(learns) #get mean learning rate
    println("learned after ", learn, " trials")
    return learn
end

function test_params(params1, params2, name1, name2; fname="test", n=100, ntrials=200)
    #test the effect of parameters on the performance of the system
    #paramsi, name1 specify the set of parameters used and which parameter to vary
    N, M = length(params1), length(params2)
    learns = zeros(N,M)
    for (i,p1) in enumerate(params1)
        for (j,p2) in enumerate(params2) #consider every parameter combination
            if (name1, name2) == ("delta", "lambda")
                learn = test_performance(n=n, ntrials=ntrials, Plot=false,
                                delta=p1, lambda=p2)
            elseif (name1, name2) == ("alpha", "beta")
                learn = test_performance(n=n, ntrials=ntrials, Plot=false,
                                alpha=p1, beta=p2)
            elseif (name1, name2) == ("mp", "l")
                learn = test_performance(n=n, ntrials=ntrials, Plot=false,
                                        mp=p1, l=p2)
            end
            learns[i,j] = learn #store performance
            print(name1,": ", p1, "  ",name2,": ",p2)
        end
    end
    heatmap(learns, params2, params1; xlab=name2, ylab=name1,
        filename="figures/"*fname*"_heat.png", Title=name1*" "*name2) #plot
    writedlm(fname*"_param1.dlm", params1)
    writedlm(fname*"_param2.dlm", params2)
    writedlm(fname*"_results.dlm", learns)
end

function compare_1_2_layer(N=500000, n = 100, ntrials = 5000, fname="compare_layers")
    #compare the performance of one- and two-layer Anderson networks
    perfs1 = zeros(ntrials, n)
    perfs2 = zeros(ntrials, n)
    for i = 1:n #average performance of n trials
        println("new i:", i)
        perf1 = run_sim(;ntrials = ntrials, mc=mc, mp=mp, l=l, muc=muc, mup=mup,
                Plot=false, alpha=1, beta=0.2, gamma=0.9, Print=false,
                layer2 = true, layer1=true) #run 1-layer simulation
        perfs1[:,i] = perf1
        perf2 = run_sim(;ntrials = ntrials, mc=mc, mp=mp, l=l, muc=muc, mup=mup,
                Plot=false, alpha=1, beta=0.2, gamma=0.9, Print=false,
                layer2 = true, layer1=false) #run 2-layer simulation
        perfs2[:,i] = perf2
    end
    #average performance
    avgperf1 = mean(perfs1, dims=2)
    avgperf2 = mean(perfs2, dims=2)

    #plot result
    figure(figsize = (6,2))
    semilogy()
    plot(1:ntrials, avgperf1, "r-")
    plot(1:ntrials, avgperf2, "b--")
    legend(["1 layer", "2 layers"])
    xlabel("trial")
    ylabel("time (s)")
    savefig("figures/"*fname*".png", bbox_inches = "tight")
    close()
end

run_sim(fname="test", N=500000, ntrials=100, all=true) #run Barto simulation
test_performance(n=100, fname="real") #average over 100 simulations

test_params([0.4;0.6;0.8;0.9;0.99],[0.4;0.6;0.8;0.9;0.99],"delta", "lambda",
            fname="test_delta_lambda_0") #test the effect of delta and lambda
test_params([0.001;1; 10; 1000; 10000],[0.05;0.25;0.5;0.75;1.0],"alpha", "beta",
            fname="test_alpha_beta") #alpha and beta
test_params([0.01;0.03;0.1;0.3;1],[0.05;0.2;0.5;2;10],"mp", "l", fname="test_mp_l")

compare_1_2_layer() #compare Anderson networks with 1 and 2 layers
