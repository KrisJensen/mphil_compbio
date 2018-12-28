using Random, PyCall, PyPlot, LinearAlgebra, DelimitedFiles
using StatsBase

function viterbi(;emitted = "emitted.txt", A="A.txt", B="B.txt",
                u0 = "u0.txt")


    emitted, A, B, u0 = get_mats(emitted, A, B, u0)

    L = length(emitted); Nh, No = size(B)
    V = zeros(L, Nh)
    traceback = zeros(L, Nh); traceback[1,:] = zeros(Nh)


    omegas = zeros(L, Nh)
    traceback = zeros(L, Nh); traceback[1,:] = zeros(Nh)
    omegas[1,:] = log10.(u0) .+ log10.(B[:, Int(emitted[1])])
    for n = 2:L
        #calculate probabilities of possible transitions
        #println(size(log10.(B[:, Int(emitted[n])])'), " ", size(omegas[n-1, :]), " ", size(log10.(A)))
        params = (log10.(A) .+ log10.(B[:, Int(emitted[n])])' ).+ omegas[n-1, :]
        #println(params)
        a, b = findmax(params, dims=1) #find optimum step to each state
        omegas[n,:] = a #store new probability
        traceback[n,:] = [b[1][1], b[2][1]] #store state from which we came
    end

    for n in 1:L println(omegas[n,:],"  ",traceback[n, :]) end

    path = zeros(L)
    # The highest probability
    opt, ind = findmax(omegas[end,:])
    path[L] = ind
    for n in reverse(1:(L-1))
        path[n] = traceback[(n+1), Int(path[n+1])]
    end

    println("Most probable path:\n",ind)
    for state in path println(state) end
    println("max prob: ", opt)

    return path, opt
end
