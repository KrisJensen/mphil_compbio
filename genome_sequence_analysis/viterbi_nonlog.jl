using Random, PyCall, PyPlot, LinearAlgebra, DelimitedFiles
using StatsBase

function viterbi(;emitted = "emitted.txt", A="A.txt", B="B.txt",
                u0 = "u0.txt")

    #implement with logs

    emitted, A, B, u0 = get_mats(emitted, A, B, u0)

    L = length(emitted); Nh, No = size(B)
    V = zeros(L, Nh)
    traceback = zeros(L, Nh); traceback[1,:] = zeros(Nh)


    V[1,:] = u0 .* B[:, Int(emitted[1])] #store viterbi probabilities
    for n in 2:L
        for i in 1:Nh
            max_tr_prob = V[(n-1), 1] * A[1, i]
            prev_st_selected = 1

            for j in 2:Nh
                tr_prob = V[(n-1),j] * A[j, i]
                if tr_prob > max_tr_prob
                    max_tr_prob = tr_prob
                    prev_st_selected = j
                end
            end

            max_prob = max_tr_prob * B[i, Int(emitted[n])]
            V[n,i] = max_prob; traceback[n,i] = prev_st_selected

        end
    end

    for n in 1:L println(V[n,:],"  ",traceback[n, :]) end
    path = zeros(L)
    # The highest probability
    opt, ind = findmax(V[end,:])
    path[L] = ind
    println("Most probable path:\n",ind)
    for n in reverse(1:(L-1))
        path[n] = traceback[(n+1), Int(path[n+1])]
        #println(path[n])
    end
    println("max prob: ",opt)

    return path, opt

end


viterbi()

#path1, opt1 = viterbi(emitted = "GC_seq.txt")

#path2, opt2 = viterbi(emitted = "GC_seq.txt", A = "Abw.txt", B="Bbw.txt", u0="ubw.txt")
