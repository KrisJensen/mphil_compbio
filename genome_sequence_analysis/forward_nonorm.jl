#implement forward, backwards and Baum-Welch algorithms

using Random, PyCall, PyPlot, LinearAlgebra, DelimitedFiles
using StatsBase

function get_mats(emitted, A, B, u0)
    #if provided as strings, returns as matrices
    if (typeof(A) == String) A = readdlm(A) end
    if (typeof(B) == String) B = readdlm(B) end
    if (typeof(u0) == String) u0 = readdlm(u0) end
    if (typeof(emitted) == String) emitted = readdlm(emitted) end
    return emitted, A, B, u0
end

function forward(;emitted = "emitted.txt", A="A.txt", B="B.txt",
                u0 = "u0.txt", Print = false)

    #load matrices if provided as files
    emitted, A, B, u0 = get_mats(emitted, A, B, u0)

    Nh, No = size(B) #number of hidden and observed states
    states_h = 1:Nh ; states_o = 1:No #state spaces
    L = size(emitted)[1]


    #fwd = zeros(length(emitted), Nh) #create array of alpha values
    fwdh = zeros(length(emitted), Nh) #create array of alpha values
    #rows are steps in the chain, columns are hidden states
    cs = zeros(length(emitted)) #store scaling factors
    #alpha_prev = zeros(Nh)
    alphah_prev = zeros(Nh) #scaled
    #initialize alpha_0(j) as u0(j)*B[j, V0]

    c1 = sum( u0 .* B[:, Int(emitted[1])] ); cs[1] = c1 #scaling factor
    for state in states_h
        #alpha_prev[Int(state)] = B[Int(state), Int(emitted[1])] * u0[Int(state)]
        alphah_prev[Int(state)] = B[Int(state), Int(emitted[1])] * u0[Int(state)] / c1
    end

    #fwd[1,:] = alpha_prev
    fwdh[1,:] = alphah_prev #store scaled values

    for (n, obs_n) in enumerate(emitted[2:end])
        #for each observed element in the chain

        #global alpha_curr = zeros(Nh) #initialize alphas
        deltas = zeros(Nh) #delta_ni = cn*alphah_ni = B[i, Vn]*sum_j[ alphah_n-1(j) * A[j, i] ]

        for state in states_h
            #update alpha_n(i) as B[i, Vn]*sum_j[ alpha_n-1(j) * A[j, i] ]
            #println(alpha_prev, " ", A[:,Int(state)], " ", alpha_prev .* A[:,Int(state)])
            #alpha_curr[state] = B[Int(state), Int(obs_n)] .*
            #                    sum(alpha_prev .* A[:,Int(state)])
            deltas[state] = B[Int(state), Int(obs_n)] .*
                                sum(alphah_prev .* A[:,Int(state)])
        end

        cn = sum(deltas) #sum_i(alphah_ni * cn) = cn * 1 since alpha_h normalized
        alphah_curr = deltas ./ cn #normalized

        #fwd[n+1,:] = alpha_curr #store alpha
        fwdh[n+1,:] = alphah_curr #store alpha hat
        #alpha_prev = alpha_curr #store as previous for next iteration
        alphah_prev = alphah_curr #store alphah
        cs[n+1] = cn #store scaling factors
        Print && println("new n: ", n+1, "  alpha: ", alpha_curr) #print result
    end

    #p = sum(alpha_curr)  #total p is sum of final alphas
    log10p = sum( log10.( cs ) ) #P = product_n(cn)
    Print && println("probability is ", 10^log10p, "  log10: ", log10p)


  # backward part of the algorithm

    #beta_prev = repeat([1], Nh)
    betah_prev = repeat([1], Nh)
    #bwd = zeros(length(emitted), Nh) #create array of alpha values
    bwdh = zeros(length(emitted), Nh) #create array of alpha values

    #bwd[1,:] = beta_prev
    bwdh[1,:] = betah_prev

    for (n, obs_Lm) in enumerate(reverse(emitted[2:end]))

        #global beta_curr = zeros(Nh)
        epsilons = zeros(Nh) #epsilon_ni = cn+1*betah_ni = sum_j[ B[j, Vn+1]* betah_n+1(j) * A[i, j] ]

        for state in states_h
            #update beta_n-1(i) as sum_j[ B[j, Vn]*beta_n(j) * A[i, j] ]
            #beta_curr[state] = sum( B[:, Int(obs_Lm)] .*
            #                    beta_prev .* A[Int(state),:] )
            epsilons[state] = sum( B[:, Int(obs_Lm)] .*
                                betah_prev .* A[Int(state),:] )
        end

        betah_curr = epsilons ./ cs[Int(L-n+1)] #normalized

        #bwd[n+1,:] = beta_curr #store beta
        bwdh[n+1,:] = betah_curr #store beta hat
        #beta_prev = beta_curr #store as previous for next iteration
        betah_prev = betah_curr #store betah
        Print && println("new n: ", L-n, "  betah: ", betah_curr) #print result
    end

    #bwd = reverse(bwd, dims=1)
    bwdh = reverse(bwdh, dims=1)

    #p_bwd = sum( bwd[1,:] .* B[:, Int(emitted[1])] .* u0 )
    #Print && println("p_bkwd: ", p_bwd)

    #P = sum(alpha_n.*beta_n) for all n!!

    return fwdh, bwdh, cs, log10p
end



#fwd_GC, fwdh_GC, cs_GC = forward(emitted = "GC_seq.txt")

function estimate_params(emitted, fwdh, bwdh, cs, A0, B0; Print = false)

    #everything is ratios of products of alphas and betas
    #so the scaling cancels (alpha_i)*(beta_i) scaled by prod(cs)

    L = size(fwdh)[1]; Nh = size(fwdh)[2]; No = size(B)[2]
    A1 = copy(A0); B1 = copy(B0)

    gammas = fwdh .* bwdh ./ sum(fwdh .* bwdh, dims = 2)

    #println(L, " ", Nh, " ", No)

    xis = zeros((L-1), Nh, Nh)
    #println(size(xis), " ", size(fwdh), " ", size(bwdh), " ", size(emitted))

    for n in  1:(L-1)
        for i in 1:Nh
            for j in 1:Nh
                #println(size(xis))
                xis[n, i, j] =
                    fwdh[n, i] * A[i,j] * bwdh[n+1, j] * B[j, Int(emitted[n+1])]

            end
        end
        xis[n, :, :] /= sum(xis[n, :, :])
        Print && println(xis[n, :, :])
    end

    u0 = gammas[1,:]
    for i in 1:Nh
        for j in 1:Nh
            A1[i, j] = sum( xis[:, i, j] ) / sum( gammas[1:(end-1),i] )
        end
        for k in 1:No
            B1[i, k] = sum( gammas[reshape(emitted .== k, L), i] ) / sum(gammas[:,i])
        end
    end

    #A1 = A1 ./ sum(A1, dims=2) #prevent drift in normalization from iterations and long chain

    return A1, B1, u0
end

function baum_welch(;emitted = "emitted.txt", A="A.txt", B="B.txt",
                u0 = "u0.txt", thresh = 0.0001, maxiter = 10000)
    #works with unscaled values; now get it working with scaled values


    #alpha_n = alphah_n * prod(cs[1:n])
    #beta_n = betah_n * prod(cs[(n+1):end])

    emitted, A_0, B_0, u0_0 = get_mats(emitted, A, B, u0)

    error = thresh+1

    niter = 0

    println(A_0)

    while (error > thresh) & (niter < maxiter)
        niter += 1
        fwdh, bwdh, cs, logp = forward(emitted=emitted, A=A_0, B=B_0, u0=u0_0)
        A_1, B_1, u0_1 = estimate_params(emitted, fwdh, bwdh, cs, A_0, B_0)
        error = norm( hcat(A_1, B_1, u0_1) - hcat(A_0, B_0, u0_0) )
        println("n: ", niter, "  error: ", error, "  logp: ", logp, "  u0 ", u0_1)
        #println(A_0, "   ", A_1)
        #println(B_0, "   ", B_1)
        #println(u0_0, "   ", u0_1)
        A_0, B_0, u0_0 = copy(A_1), copy(B_1), copy(u0_1)
    end

    println(A_0)

    fwdh, bwdh, cs, logp = forward(emitted=emitted, A=A_0, B=B_0, u0=u0_0)
    println("Final log10 Probability: ", logp)

    return A_0, B_0, u0_0, logp

end

function test_welch(N = 100)

    logps = zeros(N)
    logp_max = -Inf
    Am, Bm, um = 0, 0, 0

    for i in 1:N
        writemats(K = 2, filename = "test")
        Btest = rand(2,5) #normalize rows!
        Btest ./= sum(Btest,dims = 2)

        A, B, u0, logp = baum_welch(A = "testmat.txt", B = Btest,
                        u0 = "testvec.txt", emitted = "GC_seq.txt")

        logps[i] = logp

        if logp > logp_max
            println("New max logp: ", logp)
            Am, Bm, um = A, B, u0
            logp_max = logp
        end
    end

    return logps, Am, Bm, um
end
