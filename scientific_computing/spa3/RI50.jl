#functions for calculating regularity indices and similarities
include("plots.jl")
include("dmin2d.jl")

function calc_RI(points; Print = true)
    #given an array of points, calculates regularity index as
    # mean(dmin(i))/sd(dmin(i)) (dmin(i) is distance from i to nearest neighbor)
    n = size(points)[1] #number of points
    n < 3 && return(Inf) #sd is zero in this case

    dmins = zeros(n)
    for i = 1:n #for each point, calculate distance to all other points
        point = points[i,1:2]
        dmin = Inf
        for j = vcat(1:i-1, i+1:n) #don't need self-distance. in julia, (1:0) is empty
            d = norm(point-points[j,1:2])
            if d < dmin
                dmin = d #find nearest neighbor distance
            end
        end
        dmins[i] = dmin
    end
    RI = mean(dmins)/std(dmins) #calculate RI
    Print && println("RI is ", RI, " mean ", mean(dmins), " std ", std(dmins))
    return(RI)
end

function run_sims(;N = 1000, n=200, m=0, s=0, xlo=0, xhi=1000, ylo=0, yhi=1000)
    #given a set of parameters, generates N dmin2d point arrays
    #and calculates regularity indices. Reports 50th largest value.

    RIs = zeros(N) #initialize RI vector
    for i in 1:N
        points = dmin2d(n, m, s, xlo, xhi, ylo, yhi,
                Print=false, Plot=false)
        #get points and calculate RI
        RIs[i] = calc_RI(points, Print=false)
    end
    sort!(RIs, rev=true) #sort in place
    println("50th RI ", RIs[50], " mean ", mean(RIs),
            " median ", median(RIs), " std ", std(RIs))
    return(RIs[50]) #return 5th value (95th percentile)
end

function vary_params(;ns = [20, 60, 100, 150, 200, 300, 400, 550, 700, 850, 1000],
    ratios = [1, 3, 6, 10, 15, 20, 27.5, 35, 42.5, 50], Plot=true)
    #Investigate effect of number of points and shape of grid
    #Keep grid area constant

    results = zeros(length(ns), length(ratios)) #initialize array for RI50 values

    for (i, n) in enumerate(ns)
        for (j, ratio) in enumerate(ratios)
            #consider every combination
            print("n: ", n, " ratio: ", ratio, ":    ")
            #run simulation
            RI = run_sims(n=n, xhi = sqrt(10^6)*ratio, yhi=sqrt(10^6/ratio))
            results[i, j] = RI
        end
    end
    #plot and write results
    Plot && heatmap(results, ratios, ns, xlab="ratio", ylab="points",
            filename="RI50s.png", Title="variation of RI50")
    writedlm("vary_params.txt", results)
    return(results)
end


function calc_sim(m, s, ref; n = 100, xlo = 0, xhi = 400, ylo = 0, yhi = 400, Print=true)
    #calculate u-score

    L = size(ref)[1] #number of points
    RIs = zeros(n); us = zeros(n) #initialize arrays for storing RI50s and us
    RIs[1] = calc_RI(ref, Print=false) #calculate reference RI50

    for i in 2:n #generate n-1 reference grids
        points = dmin2d(L, m, s, xlo, xhi, ylo, yhi,
                    Print=false, Plot=false)
        RIs[i] = calc_RI(points, Print=false) #store RI50
    end

    for i in 1:n
        #u is magnitude of difference from RI to mean of other RIs
        u = abs( RIs[i] - 1/(n-1)*( sum(RIs[1:i-1])+sum(RIs[i+1:n]) ) )
        us[i] = u
    end

    #return results; u1 most informative
    Print && println("u1: ", us[1], " mean: ", mean(us), " sd: ", std(us))
    return us[1]
end

function scan_ms(;ms=0:0.5:21, ss=0:0.15:6.3, Plot=true, ref=ref, niter=1000)
    #preliminary investigation of the effect of mean and sd on u1 from ref

    results = zeros(length(ms), length(ss)) #initialize array for results
    result_mat = zeros(0,3) #store results as sequential array
    for (i, m) in enumerate(ms)
        for (j, s) in enumerate(ss)
            #consider every combination of mean and standard deviation
            print(m, " ", s, ": ")
            u = calc_sim(m, s, ref, n=niter)
            results[i,j]=u
            result_mat = [result_mat; [m, s, u]']
        end
    end

    #store and plot results
    writedlm("scan_results.txt", result_mat)
    writedlm("scan_result_mat.txt", results)
    Plot && heatmap(results, ss, ms, ylab="mean", xlab="standard deviation",
        filename="similarity.png", Title="variation of u-score")
    #report optimum parameters
    minu = findmin(results)
    println("min u is ", minu[1], " at (",
    ms[minu[2][1]], ",", ss[minu[2][2]], ")")

    return result_mat
end

function steepest_descent(; thresh=0.01, delta=0.05, rate = 1.5,
        nlim=10000, ref=ref, niter=1000, sd = "variable")
    #use a steepest descent algorithm to find the optimum set of parameters
    #use adaptive number of iterations, delta for calculating gradient
    #and step size

    params = rand(2).*[20, 6.3] #random initial parameters
    if sd != "variable" params[2] = sd end
    print("initial params ", params, " :   ")
    u = calc_sim(params[1], params[2], ref, n=10)
    results = append!(copy(params), u)'
    err = thresh+1
    n = 0

    while u > thresh && n <= nlim
        niter = Int(round(max(5000*exp(-7*u), 20))) #adaptive n
        delta = u/2 #adaptive delta for calculating gradient
        rate = u/1 #adaptive learning rate
        #calculate u1 after step in m and sd
        if sd == "variable"
            us = [calc_sim(params[1]+delta, params[2], ref, Print=false, n=niter),
                calc_sim(params[1], params[2]+delta, ref, Print=false, n=niter)]
            params += ([u, u] - us)/delta*rate #update parameters
        else
            um = calc_sim(params[1]+delta, params[2], ref, Print=false, n=niter)
            params += ([u, u] - [um, u])/delta*rate #update parameters
        end
        print("new params ", round.(params, digits=3), " n=", niter, " :   ")
        u = calc_sim(params[1], params[2], ref, n=niter) #calculate new u
        results = [results; append!(copy(params), u)'] #store result
        n += 1
    end
    writedlm("descent_results_temp.txt", results)
    return results
end

function test_sd(;n = 10, niter = 1000, ref=ref)
    #colculate the standard deviation of u at a given point
    #for a number of simulations
    us = zeros(10)
    for i in 1:10
        us[i] = calc_sim(12, 3, ref, n=niter)
    end
    sd = std(us)
    println("sd = ", sd) #result 0.011
    return sd
end

function get_opt_means(;sds = 0:0.2:5.4)
    means = zeros(length(sds))
    for (i, sd) in enumerate(sds)
        res = steepest_descent(sd = sd)
        means[i] = res[end, 1]
    end
    writedlm("opt_ms.txt", hcat(sds, means))
    return sds, means
end
