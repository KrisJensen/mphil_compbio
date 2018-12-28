using Random, PyCall, PyPlot, LinearAlgebra, DelimitedFiles, StatsBase

function writemats(;K=3, filename = "K")
    Kmat = rand(K,K) #normalize rows!
    Kmat ./= sum(Kmat,dims = 2)
    writedlm(filename*"mat.txt", Kmat)

    Kvec = rand(K)
    Kvec = Kvec ./ sum(Kvec)
    writedlm(filename*"vec.txt", Kvec)
    return Kmat, Kvec
end

#writemats()

function MC(;N=100, A = "Kmat.txt", u0 = "Kvec.txt", Plot = false)
    A = readdlm(A)
    u0 = readdlm(u0)
    S = 1:length(u0)
    chain = zeros(N)
    chain[1] = sample(S, aweights(u0))
    for i in 2:N
        chain[i] = sample(S, aweights(A[Int(chain[i-1]),:]))
    end
    if Plot
        plot(chain[1:100], "b-")
        yticks(1:size(A)[1])
        xlabel("element")
        ylabel("state")
        savefig("MC.png")
        close()
    end
    writedlm("chain1000.txt", chain)
    return chain
end

function MLE(seq, states)
    #sequence is a sequence of state indices
    #states is a set of possible state indices

    N = length(states)
    A = zeros(N, N)

    for ind in 2:length(seq)
        A[Int(seq[Int(ind-1)]), Int(seq[Int(ind)])] += 1
    end

    println("N: ", length(seq), " transitions: ", sum(A))
    A ./= sum(A, dims = 2)

    #Aij = nij / sum_i[nij]
    #where nij is #transitions from i to j
    u0 = zeros(N)
    u0[Int(seq[1])] = 1 # we only have a single data point...
    return A, u0
end

function HMC(;N=115, A = "A.txt", B = "B.txt", u0 = "u0.txt")
    hidden = MC(N=N, A=A, u0=u0)
    B = readdlm(B)
    V = 1:size(B)[2] #number of rows
    observed = zeros(N)
    for i in 1:N
        observed[i] = sample(V, aweights(B[Int(hidden[i]),:]))
    end
    return hidden, observed

end

function plot_HMC()
    h, o = HMC()
    writedlm("emitted.txt", o)
    L = length(h)
    xs = vcat(0, repeat(1:(L-1), inner=2), L)
    figure()
    psize = 1
    plot(xs, repeat(h, inner=2), "bo", MarkerSize=psize)
    plot(xs, repeat(o, inner=2), "ro", MarkerSize=psize)
    #legend(["Hidden", "Observed"])
    yticks(1:5, 1:5)
    xlabel("state")
    ylabel("element")
    savefig("ho.png")
    #show()
    close("all")
end
