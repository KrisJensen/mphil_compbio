using Random, PyCall, PyPlot, LinearAlgebra, DelimitedFiles, StatsBase

function get_brackets(;file = "S.cerevisiae.III.fa")
    #calculate GC content of file
    #convert to number 1:5 by binning for emission
    #character count is multiple of 60
    f = open(file)
    lines = readlines(f)
    genome = ""
    for line in lines
        genome = genome*line #construct single string with genome
    end
    println(length(genome), " bases" )
    n = 1
    brackets = []
    while n <= length(genome)
        brackets = vcat(brackets,
                        genome[n:min((n+99), length(genome))])
        n += 100 #store in windows of 100 bp
    end
    close(file)
    return brackets
end

function calc_GC(brackets)
    #calculate GC content of each bin and return sequence of GC contents
    seq = zeros(length(brackets))
    for i in 1:length(seq)
        seq[i] = length( collect( #search for any of G, C, g, c
                eachmatch(r"[GCgc]", brackets[i], overlap=false) ) ) /
                length(brackets[i]) #normalize by window length
    end
    return seq
end

function anal_GC(;N=100000)
#for binning purposes, we aim to match the probability
#of getting a particular observable in the model with the probability
#of getting a particular GC content in the chromosome
#i.e. match cumulative distributions. Can also do uniform binning
 h, o = HMC(N=N)
 cum = 0
 for i in 1:5
     cum += sum(o .== i)/N
     println(i, ": ", sum(o .== i)/N, "  ", cum )
     #print cumulative distribution of HMC
 end
 bins = [0.275, 0.345, 0.405, 0.495, 1.0]
 #bins_eq = [0.325, 0.365, 0.405, 0.445, 1.0]
 for i in 1:5
     println(sum( seq .<= bins[i] )/length(seq))
     #print cumulative distribution of GC content
 end
end

function GC_to_state(seq, bins = [0.275, 0.345, 0.405, 0.495, 1.0])
    #convert from GC content to state 1:5
    #for use with HMM
    N = length(seq)
    state_seq = zeros(N) #initialize state sequence
    for i in 1:N
        for j = 1:5
            if seq[i] <= bins[j]
                state_seq[i] = j
                break
            end
        end
    end
    return state_seq
end

function plot_GCs(;viterbi = "viterbi_path_bw.txt", GC_obs = "GC_seq.txt",
                    GC_percent = "GC_percent.txt")
    #plot GC content for a subset of yeast chromosome III
    #together with hidden states inferred with the Viterbi algorithm

    Range = [100, 300]
    viterbi = readdlm(viterbi) #read hidden states
    GC_obs = readdlm(GC_obs) #read observed GC states
    GC_percent = readdlm(GC_percent) #read GC content
    min, max = minimum(GC_percent), maximum(GC_percent)
    GC_percent = (GC_percent .- min) * 3/(max-min) #rescale range to [0:3]

    figure()
    psize = 1
    plot(viterbi, "ko", MarkerSize = psize)
    #plot(GC_obs, "bo", MarkerSize = psize)
    plot(GC_percent, "ro", MarkerSize = psize)
    xlim(Range)
    xlabel("bracket")
    ylabel("state / scaled GC content")
    legend(["Viterbi path", "GC content"])
    savefig("viterbi_path.png")
    close()
end
