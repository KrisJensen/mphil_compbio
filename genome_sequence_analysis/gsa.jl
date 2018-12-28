###code for GSA

using Random, PyCall, PyPlot, LinearAlgebra, DelimitedFiles, StatsBase
include("forward.jl")
include("viterbi.jl")
include("markov_chains.jl")
include("anal_GC.jl")

    #1.
#Kmat, Kvec = writemats(K = 3)
#chain = MC( N = 1000, Plot=true ) #generate example chain

    #2.
#A, u0 = MLE(chain, 1:3) #infer parameters (have to repeat this to get u0)

    #3.
#construct model from assignment
A = [0.8 0.2; 0.1 0.9]; writedlm("A.txt", A)
u0 = [0.5 ; 0.5]; writedlm("u0.txt", u0)
B = [0.2 0.5 0.2 0.1 0; 0 0.1 0.4 0.4 0.1]; writedlm("B.txt", B)
#plot_HMC() #generate and plot chain for hidden markov model

    #4.
#use forward algorithm to get probability of emitted sequence
#fwdh, bwdh, cs, log10p = forward(emitted = "emitted.txt"))

    #5.
bs = get_brackets() #load yeast genome and divide into 100 bp windows
seq = calc_GC(bs) #calculate GC content of bins

#plot cumulative distribution
#PyPlot.plt[:hist](seq, bins = 0.05:0.01:0.95, cumulative = true, density=true)
    #xlabel("GC content")
    #ylabel("Frequency")
    #savefig("cumulative_GC.png")
#close()

#state_seq = GC_to_state(seq) #infer discrete states 1:5
#writedlm("GC_seq.txt", state_seq); writedlm("GC_percent,txt", seq)
#get log likelihood of GC sequence
#fwdh, bwdh, cs, log10p = forward(emitted = "GC_seq.txt"))

    #6.
#Use Baum-Welch to infer parameters with initial parameters from assignment
#A, B, u0, logp = baum_welch(emitted = "GC_seq.txt")
#use uniform binning
#Aunif, Bunif, u0unif, logpunif = baum_welch(emitted = "GC_seq_unif.txt")

#test 100 different initial parameters and find best possible parameters
#logps, Am, Bm, um = test_welch()
#use uniform binning
#logps, Am, Bm, um = test_welch(emitted = "GC_seq_unif.txt")

    #7.
#p, o = viterbi() #test viterbi for N=115 simulated HMM
#generate Viterbi path for Baum-Welch-inferred parameters
#path, opt = viterbi(emitted = "GC_seq.txt", A = "Abw.txt", B="Bbw.txt", u0="ubw.txt")
#writedlm("viterbi_path_bw.txt", path)

#plot_GCs() #plot viterbi path together with GC content
