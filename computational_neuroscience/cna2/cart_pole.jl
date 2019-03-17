##code for the art pole balancing problem

using PyCall, PyPlot, LinearAlgebra, Distributions,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats


###TRY HAVING POSITIVE GRAVITATIONAL CONSTANT!

#Physical parameters (vary e.g. mass?)
g = -9.8 #m/s2, acceleration due to gravity
g = 9.8
mc = 1.0 #kg, mass of cart
mp = 0.1 #kg, mass of pole,
l = 0.5 #m, half-pole length,
muc = 0.0005#, coefficient of friction of cart on track,
mup = 0.000002#, coefficient of friction of pole on cart,
Ft = 10.0 # p/m 10 newtons, force applied to cart's center of mass at time t.

#model parameters (vary some of these)
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

    if all

    ts = (1:N)*dt
    figure(figsize = (5,3))
    plot(ts, xs)
    plot(ts, repeat([-2.4], N), "k--")
    plot(ts, repeat([2.4], N), "k--")
    xlabel("time (s)")
    ylabel("x")
    savefig("figures/"*fname*"_xs.png", bbox_inches = "tight")
    close()

    figure(figsize = (5,3))
    plot(ts, thetas)
    plot(ts, repeat([-12], N), "k--")
    plot(ts, repeat([12], N), "k--")
    xlabel("time (s)")
    ylabel("theta")
    savefig("figures/"*fname*"_thetas.png", bbox_inches = "tight")
    close()

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

    newvs = zeros(3, 6)
    for (i,thetadot) = enumerate([-100;0;100])
        for (j,theta) = enumerate([-9; -3.5; -0.5; 0.5; 3.5; 9])
            temps = zeros()
            for x = [-1.6;0;1.6]
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

    ntrial = length(performance)
    figure(figsize = (5,3))
    plot(1:ntrial, performance)
    xlabel("trial")
    ylabel("time (s)")
    savefig("figures/"*fname*"_performance.png", bbox_inches = "tight")
    close()

end

function smooth_performance(performance, dt)
    N = length(performance)
    newperf = zeros(N)
    for i = 1:4
        newperf[i] = mean(performance[1:i])
    end
    for i = 5:N
        newperf[i] = mean(performance[i-4:i])
    end
    return newperf*dt #real time
end

function get_state(state)
    x, theta, xdot, thetadot = state

    k1 = sum(x .> [-0.8, 0.8])
    k2 = sum(theta .> [-6,-1,0,1,6])
    k3 = sum(xdot .> [-0.5, 0.5])
    k4 = sum(thetadot .> [-50, 50])
    k = 54*k1 + 9*k2 + 3*k3 + k4 + 1
    return k #index of non-zero element

    vec = zeros(162)
    vec[k] = 1
    return vec
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
    #println(alpha, " ", r, " ", es)
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
    x, theta, xdot, thetadot = state
    theta, thetadot = 2*pi/360*theta, 2*pi/360*thetadot
    Ft = 10*yt
    #work in radians
    dthetadot = ( g*sin(theta) +
                cos(theta)*( ( -Ft-mp*l*thetadot^2*sin(theta) + muc*sign(xdot) )/( mc + mp ) ) -
                mup*thetadot/(mp*l) ) / (
                l*( 4/3 - ( mp* (cos(theta)^2) )/( mc+mp ) ))
    dxdot = (  Ft + mp*l*(thetadot^2*sin(theta) - dthetadot*cos(theta)) -
            muc*sign(xdot) ) / ( mc+mp )
    thetadot += dthetadot*dt
    xdot += dxdot*dt
    theta += thetadot*dt
    x += xdot*dt
    return [x; 360/(2*pi)*theta; xdot; 360/(2*pi)*thetadot]
end

function get_zs(k, D)
    D[i,k]

end

function update_2layer(state, ws, vs, fs, D, cs, A, ptm, ktm; d=dnorm,
                        alpha=alpha, lambda=lambda, gamma=gamma, delta=delta, beta=beta,
                        l=l, mc=mc, mp=mp, muc=muc, mup=mup)

    n = 162
    #ASE
    k = get_state(state)
    zs = 1 ./ (1 .+ exp.(-D[:,k])) #vector
    P = ws[k] + fs' * zs
    P = 1 / (1 + exp(-P)) #probability of right
    #println(P)
    if rand() < P
        yt = 1
    else
        yt = -1
    end
    #ACE
    qs = 1 ./ (1 .+ exp.(-A[:,k])) #vector
    #if ktm == 0 qst1t1 = zeros(n) else qst1t1 = 1 ./ (1 .+ exp.(-A[:,ktm])) end
    pt = vs[k] + cs' * qs
    #if ktm == 0 pt1t1 = 0 else pt1t1 = vs[ktm] + cs' * qst1t1 end

    if state == [0;0;0;0]
        rhatt = 0
    else
        rhatt = gamma*pt - ptm #calculate expected reward at time t
        #rhatt = gamma*pt - pt1t1
    end

    ws[k] += alpha*rhatt*(max(yt,0) - P)
    fs += alpha * rhatt * (max(yt,0) - P)*zs
    D[:,k] += ((0.2*alpha) * rhatt * zs) .* (1 .- zs) .* sign.(fs) * (max(yt,0) - P)

    vs[k] += beta*rhatt
    A[:,k] += ((0.25*beta)*rhatt*qs) .* (1 .- qs) .* sign.(cs)
    #A[:,k] += ((0.25*beta)*rhatt*qst1t1) .* (1 .- qst1t1) .* sign.(cs)
    cs += beta * rhatt * qs
    #cs += beta * rhatt * qst1t1


    state = update_cart(state, yt, l=l, mc=mc, mp=mp, muc=muc, mup=mup) #update states to t+1
    return state, ws, vs, fs, D, cs, A, pt, k, yt, P, zs, qs
end

function update_system(state, ws, vs, es, xbars, ptm; d=dnorm,
                        alpha=alpha, lambda=lambda, gamma=gamma, delta=delta, beta=beta,
                        l=l, mc=mc, mp=mp, muc=muc, mup=mup)
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

function run_sim(;fname="test", N = 500000, ntrials = 100, mc=mc, mp=mp, l=l, muc=muc, mup=mup, layer2 = false,
                alpha=alpha, beta=beta, gamma=gamma, lambda=lambda, delta=delta, d=dnorm, Plot=true, Print=true,
                all = false)
    n = 162
    pt, r, k = 0, 0, 0 #no reward or predicted reward
    ws, vs, es, xbars, fs, cs = zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n)
    A, D = zeros(n,n), zeros(n,n)
    state = [0;0;0;0]
    ntrial, lasttrial, performance = 0, 0, zeros(ntrials)
    thetas = zeros(N)
    xs = zeros(N)
    Print && println("\nnew simulation")
    for i = 1:N
        if !layer2
        state, ws, vs, es, xbars, pt, k = update_system(state, ws, vs, es, xbars, pt,
                                                    l=l, mc=mc, mp=mp, muc=muc, mup=mup,
                                                    alpha=alpha, beta=beta, gamma=gamma,
                                                    lambda=lambda, delta=delta)
        else
        state, ws, vs, fs, D, cs, A, pt, k, yt, P, zs, qs = update_2layer(state, ws, vs, fs, D, cs, A, pt, k,
                                                        l=l, mc=mc, mp=mp, muc=muc, mup=mup,
                                                        alpha=alpha, beta=beta, gamma=gamma,
                                                        lambda=lambda, delta=delta)
        end
        xs[i] = state[1]
        thetas[i] = state[2]

        if abs(state[1]) >= 2.4 || abs(state[2]) >= 12 #we have failed
            ntrial += 1
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
            #vs = vs + beta* rhatt *xbars #learn from real error signal
            #xbars *= lambda #update xbars
            #es = zeros(n) #reset eligibility trace

            if !layer2
                ws = update_ws(ws, rhatt, es, alpha) #update ASE weights t+1
                vs = update_vs(vs, rhatt, xbars, beta) #update ACE weights t+1
                es = zeros(n) #update_es(es, k, yt, delta) #update ASE eligibility traces t+1
                xbars = update_xbars(xbars, k, lambda)
            else
                ws[k] += alpha*rhatt*(max(yt,0) - P)
                fs += alpha * rhatt * (max(yt,0) - P)*zs
                D[:,k] += ((0.2*alpha)*rhatt*zs) .* (1 .- zs) .* sign.(fs) * (max(yt,0) - P)
                vs[k] += beta*rhatt
                cs += beta * rhatt * qs
                A[:,k] += ((0.25*beta)*rhatt*qs) .* (1 .- qs) .* sign.(cs)
            end

            #update 1 step
            state, r = [0;0;0;0], 0#-1 #reset system and set r = -1 (negative reward)
            performance[ntrial] = i-lasttrial
            lasttrial = i
            if ntrial == ntrials break end
        else
            r = 0 #no reinforcement
        end
    end

    if ntrial < ntrials
        performance[ntrial+1:end] .= max(performance[ntrial], N-lasttrial)
    end
    performance = smooth_performance(performance, dt)
    Plot && plot_poles(N, xs, thetas, performance, dt, fname, vs; all=all)
    println(vs)
    return performance
end

function test_performance(;n=10, ntrials=100, mc=mc, mp=mp, l=l, muc=muc, mup=mup, fname="test", Plot=true,
                alpha=alpha, beta=beta, gamma=gamma, lambda=lambda, delta=delta, d=dnorm)
    perfs = zeros(ntrials, n)
    learns = zeros(n)
    for i = 1:n
        perf = run_sim(;ntrials = ntrials, mc=mc, mp=mp, l=l, muc=muc, mup=mup,Plot=false,
                alpha=alpha, beta=beta, gamma=gamma, lambda=lambda, delta=delta, Print=true)
        perfs[:,i] = perf
        learn = sum(perf .< 500000/10*0.02)
        learns[i] = learn
    end

    avgperf = mean(perfs, dims=2)

    if Plot
        figure(figsize = (5,3))
        plot(1:ntrials, avgperf)
        xlabel("trial")
        ylabel("time (s)")
        savefig("figures/"*fname*"_meanperformance.png", bbox_inches = "tight")
        close()
    end
    learn = sum(avgperf .< 500000/10*0.02)
    #learn = avgperf[end]
    learn = mean(learns)
    println("learned after ", learn, " trials")
    return learn
end

function test_params(params1, params2, name1, name2; fname="test", n=50, ntrials=200)
    N, M = length(params1), length(params2)
    learns = zeros(N,M)
    for (i,p1) in enumerate(params1)
        for (j,p2) in enumerate(params2)
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
            learns[i,j] = learn
            print(name1,": ", p1, "  ",name2,": ",p2)
        end
    end
    heatmap(learns, params2, params1; xlab=name2, ylab=name1,
        filename="figures/"*fname*"_heat.png", Title=name1*" "*name2)
    writedlm(fname*"_param1.dlm", params1)
    writedlm(fname*"_param2.dlm", params2)
    writedlm(fname*"_results.dlm", learns)
end

run_sim(fname="learn2layer", N=101, ntrials=100, all=true, layer2=true, alpha=1, beta=0.2, gamma=0.9)
#run_sim(fname="test", N=500000, ntrials=100, all=true)
#run_sim(fname="l05", l = 0.5, N=500000, ntrials=500, layer2 = true)
#run_sim(fname="l01", l = 0.1, N=50000)
#test_performance(n=100, fname="real")

#test_params([0.4;0.6;0.8;0.9;0.99],[0.4;0.6;0.8;0.9;0.99],"delta", "lambda", fname="test_delta_lambda_0")
#test_params([0.1;0.7; 10; 1000; 10000],[0.0;0.25;0.5;0.75;1.0],"alpha", "beta", fname="test_alpha_beta")
#test_params([0.01;0.03;0.1;0.3;1],[0.05;0.2;0.5;2;10],"mp", "l", fname="test_mp_l")
