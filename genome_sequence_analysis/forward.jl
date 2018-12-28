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

    fwdh = zeros(length(emitted), Nh) #create array of alpha values
    #rows are steps in the chain, columns are hidden states
    cs = zeros(length(emitted)) #store scaling factors
    alphah_prev = zeros(Nh) #scaled
    #initialize alpha_0(j) as u0(j)*B[j, V0]

    c1 = sum( u0 .* B[:, Int(emitted[1])] ); cs[1] = c1 #scaling factor
    for state in states_h
        alphah_prev[Int(state)] = B[Int(state), Int(emitted[1])] * u0[Int(state)] / c1
    end
    fwdh[1,:] = alphah_prev #store scaled values

    for (n, obs_n) in enumerate(emitted[2:end])
        #for each observed element in the chain
        deltas = zeros(Nh) #delta_ni = cn*alphah_ni = B[i, Vn]*sum_j[ alphah_n-1(j) * A[j, i] ]

        for state in states_h
            deltas[state] = B[Int(state), Int(obs_n)] .*
                                sum(alphah_prev .* A[:,Int(state)])
        end

        cn = sum(deltas) #sum_i(alphah_ni * cn) = cn * 1 since alpha_h normalized
        alphah_curr = deltas ./ cn #normalized

        fwdh[n+1,:] = alphah_curr #store alpha hat
        alphah_prev = alphah_curr #store alphah
        cs[n+1] = cn #store scaling factors
        Print && println("new n: ", n+1, "  alpha: ", alphah_curr) #print result
    end

    #p = sum(alpha_curr)  #total p is sum of final alphas
    log10p = sum( log10.( cs ) ) #P = product_n(cn)
    Print && println("probability is ", 10^log10p, "  log10: ", log10p)


  # backward part of the algorithm

    betah_prev = repeat([1], Nh)
    bwdh = zeros(length(emitted), Nh) #create array of alpha values
    bwdh[1,:] = betah_prev

    for (n, obs_Lm) in enumerate(reverse(emitted[2:end]))

        epsilons = zeros(Nh) #epsilon_ni = cn+1*betah_ni = sum_j[ B[j, Vn+1]* betah_n+1(j) * A[i, j] ]

        for state in states_h
            epsilons[state] = sum( B[:, Int(obs_Lm)] .*
                                betah_prev .* A[Int(state),:] )
        end

        betah_curr = epsilons ./ cs[Int(L-n+1)] #normalized

        bwdh[n+1,:] = betah_curr #store beta hat
        betah_prev = betah_curr #store previous beta hat
        Print && println("new n: ", L-n, "  betah: ", betah_curr) #print result
    end

    bwdh = reverse(bwdh, dims=1)
    Print && println( sum(fwdh .* bwdh, dims=2) ) #P/P = 1 for all n

    return fwdh, bwdh, cs, log10p
end


function estimate_params(emitted, fwdh, bwdh, cs, A0, B0; Print = false)
    #given alphas, betas and emitted sequence, calculate expected A, B, u0
    #everything is ratios of products of alphas and betas
    #so the scaling cancels
    #(alpha_i)*(beta_i) scaled by prod(cs)

    L = size(fwdh)[1]; Nh = size(fwdh)[2]; No = size(B0)[2]
    A1 = copy(A0); B1 = copy(B0)

    #gamma_ni = alpha_n(i)*beta_n(i) / sum_i(alpha_n(i) * beta_n(i))
    gammas = fwdh .* bwdh ./ sum(fwdh .* bwdh, dims = 2)

    xis = zeros((L-1), Nh, Nh) #initialize array

    for n in  1:(L-1)
        #transitions from n to n+1
        for i in 1:Nh
            for j in 1:Nh
                xis[n, i, j] =
                    fwdh[n, i] * A0[i,j] * bwdh[n+1, j] * B0[j, Int(emitted[n+1])]

            end
        end
        xis[n, :, :] /= sum(xis[n, :, :]) #normalize xis
        Print && println(xis[n, :, :])
    end

    u0 = gammas[1,:] #get initial distribution
    for i in 1:Nh
        for j in 1:Nh
            #update transition matrix as P(xn=i, xn+1=j) / P(xn=i)
            A1[i, j] = sum( xis[:, i, j] ) / sum( gammas[1:(end-1),i] )
        end
        for k in 1:No
            #update emission matrix as p(xn=i, bn=k) / p(xn=i)
            B1[i, k] = sum( gammas[reshape(emitted .== k, L), i] ) / sum(gammas[:,i])
        end
    end

    return A1, B1, u0
end

function baum_welch(;emitted = "emitted.txt", A="A.txt", B="B.txt",
                u0 = "u0.txt", thresh = 0.0001, maxiter = 10000)
    #given an emitted sequence and initial parameters
    #calculates optimum parameters using the Baum Welch algorithm
    #can also provide convergence threshold and max number of iterations

    #load matrices if provided as strings
    emitted, A_0, B_0, u0_0 = get_mats(emitted, A, B, u0)
    error = thresh+1
    niter = 0

    println(A_0)

    while (error > thresh) & (niter < maxiter)
        niter += 1
        #calculate alphas and betas from forward algorithm
        fwdh, bwdh, cs, logp = forward(emitted=emitted, A=A_0, B=B_0, u0=u0_0)
        #estimate expected parameters from emitted sequence and alphas+betas
        A_1, B_1, u0_1 = estimate_params(emitted, fwdh, bwdh, cs, A_0, B_0)
        #check for convergence
        error = norm( hcat(A_1, B_1, u0_1) - hcat(A_0, B_0, u0_0) )
        println("n: ", niter, "  error: ", error, "  logp: ", logp, "  u0 ", u0_1)
        A_0, B_0, u0_0 = copy(A_1), copy(B_1), copy(u0_1)
    end

    println(A_0)

    fwdh, bwdh, cs, logp = forward(emitted=emitted, A=A_0, B=B_0, u0=u0_0)
    #calculate converged parameters
    println("Final log10 Probability: ", logp)
    return A_0, B_0, u0_0, logp
end

function test_welch(;N = 100, emitted = "GC_seq.txt")
    #generate N sets of random initial parameters
    #and run Baum-Welch algorithm with these

    logps = zeros(N) #initialize some stuff
    logp_max = -Inf
    Am, Bm, um = 0, 0, 0

    for i in 1:N
        #generate initial parameters
        writemats(K = 2, filename = "test") #generates A and u0
        Btest = rand(2,5) #generate emission matrix
        Btest ./= sum(Btest,dims = 2) #normalize rows
        A, B, u0, logp = baum_welch(A = "testmat.txt", B = Btest,
                        u0 = "testvec.txt", emitted = emitted)
        logps[i] = logp #store logp
        if logp > logp_max
            #if best logp so far, store optimum parameters
            println("New max logp: ", logp)
            Am, Bm, um = A, B, u0
            logp_max = logp
        end
    end
    #return all p values and optimum parameters for best p value
    return logps, Am, Bm, um
end
