using PyCall, PyPlot, LinearAlgebra, Random, DelimitedFiles


function get_ints(spikes)
    intervals = []
    rates = []
    if spikes[1][1] < spikes[2][1]
        j = 1
    else
        j = 2
    end
    for i in 1:(length(spikes[j])-1)
        intervals = [intervals; spikes[3-j][i]-spikes[j][i]]
        rates = [rates; 1/(spikes[j][i+1]-spikes[j][i])*1000]
    end

    return [spikes[j][1:length(intervals)], intervals, rates]
end

function plotIF(times, V, spikes; interval=true, filename="test")

    figure()
    plot(times, V[:,1], "b-")
    plot(times, V[:,2], "r-")
    types = ["b-", "r-"]
    for j in 1:2
        for s in spikes[j]
            plot([s,s], [-80,0], types[j])
        end
    end
    xlim(times[1], times[end])
    xlabel("t (ms)")
    ylabel("V (mV)")
    savefig(filename*".png")
    show()
    close()

    if interval
        times, intervals, rates = get_ints(spikes)
        figure()
        plot(times, intervals)
        savefig(filename*"_int.png")
        close()

        return [times, intervals]
    end
end

function IF(V1, V2; Es=0, T=2000, dt=0.005, Pmax=0.5, ts=10, tm=20, rmgs=0.15,
            RmIet=18, EL=-70, Vt=-54, Vreset=-80, plotrange="all", filename="test")

    N = Int(round(T / dt) + 1)
    times = 0:dt:T
    V = zeros(N,2)
    Ps = zeros(N,2)
    z = zeros(N, 2)
    V[1,:] = [V1; V2]
    spikes = [ [], [] ]
    for i in 2:N
        for j in 1:2
            z[i,j] = z[i-1,j] - dt/ts * z[i-1,j]
            Ps[i,j] = Ps[i-1,j] + dt/ts * (exp(1) * Pmax*z[i,j]-Ps[i-1,j])
            V[i,j] = V[i-1, j] + dt/tm *
                (EL - V[i-1,j] -rmgs*Ps[i-1,j]*(V[i-1, j]-Es)+RmIe)
        end
        for j in 1:2
            if V[i,j] >= Vt
                z[i, 3-j] = 1
                V[i,j] = Vreset
                spikes[j] = [spikes[j];times[i]]
            end
        end
    end

    if plotrange == "none"
        a = 2
    elseif plotrange == "all"
        plotIF(times, V, spikes, filename=filename)
    else
        plotrange = Int(plotrange[1]/dt):1:Int(plotrange[2]/dt)
        plotIF(times[plotrange], V[plotrange,:], spikes, filename=filename)
    end

    return times, V, spikes, z, Ps

end

function test_intervals(;ts=10, Es=-80, rmgs = 0.15)
    ints = [[], [], [], [], [], []]
    if Es == -80
        Vis = [-71; -72; -72.35; -72.4; -73; -74]
    else
        Vis = [-70.5; -72; -74; -76; -78; -80]
    end
    for i in 1:length(Vis)
        times, V, spikes, z, Ps = IF(-70, Vis[i], Es=Es, T=3000,
            plotrange="none", ts=ts, rmgs=rmgs)
        ints[i] = plotIF(times, V, spikes)
    end
    figure()
    cols = ["b--", "r--", "g--", "y--", "c--", "m--"]
    for i in 1:length(Vis)
        plot(ints[i][1], ints[i][2])
    end
    xlabel("t (ms)")
    ylabel("interspike interval (ms)")
    legend([string(Vi) for Vi in Vis])
    savefig("intervals_ts_"*string(ts)*"_rmgs_"*string(rmgs)*"e.png")
    show()
    close()
end

function random_init(;N=75)#Es, ts, rmgs; N=5)

    Ess = [0, -80]
    tss = [5, 10, 15]
    rmgss = [0.05, 0.15, 0.30]

    n = 0
    labs = []
    ncats = length(Ess)*length(tss)*length(rmgss)-1
    ints = zeros(N, ncats)
    rates = zeros(N, ncats)
    xs = zeros(N, ncats)
    for Es in Ess
        for ts in tss
            for rmgs in rmgss
                if !(Es==-80 && ts == 15 && rmgs == 0.30) #no spikes
                    println("new params ", Es, " ", ts, " ", rmgs)
                    labs = [labs; "("*string(Es)*", "*string(ts)*", "*string(rmgs)*")"]
                    n += 1
                    #ints = []
                    #rates = []
                    for i in 1:N
                        println("new i", i)
                        v1, v2 = rand(2)*25-[80,80]
                        times, V, spikes, z, Ps = IF(v1, v2, Es=Es, T=10000,
                            plotrange="none", ts=ts, rmgs=rmgs)
                        res = get_ints(spikes)
                        is = res[2]; rs = res[3]
                        #rates = [rates; rs[end]]
                        rates[i,n] = rs[end]
                        if is[end] <= (1/rs[end])*1000/2
                            ints[i,n] = is[end]
                        else
                            ints[i,n] = (1/rs[end])*1000 - is[end]
                        end
                    end
                    xs[:,n] = rand(N)*0.8 + repeat([n-0.40], N)
                    #plot(xs, ints, "bx", MarkerSize = 1)
                end
            end
        end
    end

    fignames = ["ints_scatter", "rates_scatter"]
    ylabs = ["interspike interval (ms)", "firing rate (hz)"]
    for i in 1:2
        data = [ints, rates][i]
        figure(figsize = [11,3])
        if i == 1
            plot([0.5, ncats+0.5], [0,0], linestyle="--",color="0.1", linewidth = 0.5)
        end
        for n in 1:ncats
            plot(xs[:,n], data[:,n], "bx", MarkerSize = 1)
            plot([n+0.5, n+0.5], [-10, 100], "k--", linewidth = 0.5 )
        end
        if i == 2
            xticks(1:length(labs), labs, rotation=45, ha="right")
        else
            xticks([])
        end
        ylim([-1, maximum(data)+1])
        xlim(0.5, ncats+0.5)
        ylabel(ylabs[i])

        savefig(fignames[i]*".png", bbox_inches="tight", dpi=360)
        show()
        close()
    end
    #writedlm("ints.dlm", ints)
    #writedlm("rates.dlm", rates)
    #writedlm("xs.dlm", xs)
end

#random_init()#-80, 10, 0.15)
test_intervals(Es=0, rmgs=0.5, ts=10)
#times, V, spikes, z, Ps = IF(-65, -75,
#    Es=0, T=3001, rmgs=0.05, ts=5, plotrange=[2500,3000], filename="test")



# #potentials in mV
# EL = -70 #equlibrium leak potential
# Vt = -54 #spike threshold
# Vreset = -80 #refractory potential
#
# #times in ms
# tm = 20 #membrane time constant
# rmgs = 0.15 #synaptic conductivity times membrane resistance. units cancel
# RmIe = 18 #mV. external input current
#
# T = 500
# dt = 0.005
#ts = 10
#Pmax = 0.5

#dPs = dt / ts * ( exp(1) * Pmax*z-Ps )
#dz = - dt / ts * z

#Es = 0
