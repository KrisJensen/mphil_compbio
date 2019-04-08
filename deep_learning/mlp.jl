using PyCall, PyPlot, LinearAlgebra, Distributions,
    Random, DelimitedFiles, LaTeXStrings, StatsBase, Statistics, MultivariateStats
@pyimport pickle

##pickle functions courtesy of Rob Blackwell
#https://gist.github.com/RobBlackwell/10a1aeabeb85bbf1a17cc334e5e60acf
function mypickle(obj, filename)
    out = open(filename,"w")
    pickle.dump(obj, out)
    close(out)
 end

function unpickle(filename)
    r = nothing
    @pywith pybuiltin("open")(filename,"rb") as f begin
        r = pickle.load(f)
    end
    return r
end

function heatmap(results, xs, ys; xlab="", ylab="",
    filename="default", Title="default", vmin="default", vmax="default")
    #given a matrix of z values and lists of x,y values
    #plots a heatmap
    figure()
    if vmin == "default"
        imshow(results, cmap=ColorMap("gray_r"), vmin = minimum(results),
                vmax = maximum(results))
    else
        imshow(results, cmap=ColorMap("gray_r"), vmin = vmin,
                vmax = vmax)
    end
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

function calc_vs_zs(us; w1s=w1s, w2s=w2s, bias = false, dropout = [0 0],
                    transfer="tanh")
    #for a given input, calculates hidden and output activities
    if transfer == "tanh"
        vs = tanh.(us * w1s) #calculate hidden layer values with tanh transfer
    elseif transfer == "sigmoid"
        vs = 1 ./ (1 .+ exp.( -us * w1s )) #sigmoidal
    else println("transfer function not recognized")
    end

    if bias vs = [vs 1] end
    if dropout[2] > 0 #dropout on layer 2
        n = length(vs)
        inds = randsubseq(1:n, dropout[2])
        vs[1, inds] .= 0
        vs *= 1/(1-dropout[2]) #rescale
    end
    if transfer == "tanh"
        zs = tanh.(vs * w2s) #calculate output values with tanh
    elseif transfer == "sigmoid"
        zs = 1 ./ (1 .+ exp.( -vs * w2s )) #sigmoidal
    end
    return vs, zs
end

function backpropagate(us, vs, zs, ts; eta = 0.1, w1s=w1s, w2s=w2s, N=N, transfer="tanh")
    #update weights for second layer
    Nu = length(us)
    Nuout = size(w1s)[2] #different from Nv if bias units
    Nv = length(vs)
    Nz = length(zs)
    deltaks = zeros(Nz) #store deltas for last layer
    dw2s = zeros(Nv, Nz)
    for k = 1:Nz #for each output unit
        if transfer == "tanh"
            delta = (1 - zs[k]^2) * (ts[k]-zs[k]) #calculate delta tanh
        elseif transfer == "sigmoid"
            delta = (zs[k] - zs[k]^2) * (ts[k]-zs[k]) #sigmoidal
        end
        deltaks[k] = delta
        for j = 1:Nv #for each hidden unit
            dw2s[j, k] = eta * delta * vs[j] #update rule
        end
    end
    w2s += dw2s #modify weights

    #update weights for first layer
    deltajs = zeros(Nuout)
    dw1s = zeros(Nu, Nuout)
    for j = 1:Nuout #for each hidden unit
        if transfer == "tanh"
            delta = (1-vs[j]^2) * sum(deltaks .* w2s[j,:]) #calculate delta tanh
        elseif transfer == "sigmoid"
            delta = (vs[j]-vs[j]^2) * sum(deltaks .* w2s[j,:]) #sigmoidal
        end
        deltajs[j] = delta #store deltas
        for i = 1:Nu #for each input unit
            dw1s[i,j] = eta * delta * us[i] #update rule
        end
    end
    w1s += dw1s #modiy weights
    return w1s, w2s
end

function train_round(;ninput = ninput, train_ts=train_ts, train_us=train_us,
                w1s=w1s, w2s=w2s, eta=0.1, N=N, bias=false, dropout=[0 0], transfer="tanh")
    #trains a single epoch of the neural network using backpropagation
    nus = size(train_us)[2]
    nzs = size(train_ts)[2]
    for n = 1:ninput #for each image
        ts = reshape(train_ts[n,:], 1, nzs) #target, hvec
        us = reshape(train_us[n,:], 1, nus) #input, hvec
        if bias us = [us 1] end
        if dropout[1] > 0 #dropout on layer 1
            n = length(us)
            inds = randsubseq(1:n, dropout[1])
            us[1, inds] .= 0
            us *= 1/(1-dropout[1]) #rescale
        end
        vs, zs = calc_vs_zs(us, w1s=w1s, w2s=w2s, bias=bias, dropout = dropout, transfer=transfer) #get output
        w1s, w2s = backpropagate(us, vs, zs, ts, w1s=w1s, w2s=w2s, eta=eta, N=N, transfer=transfer) #update weights
    end
    return w1s, w2s
end

function test_all(;test_ts = test_ts, test_us = test_us,
                train_ts=train_ts, train_us = train_us,
                ntest=ntest, w1s=w1s, w2s=w2s,
                bias = false, test=true, loss=false, transfer="tanh",
                Plot = false, filename="test")
    #test performance on all test data
    #if loss, calculates a mean square error instead
    #if Plot, plots performance for individual digits
    if test #consider test data
        testinput = test_us
        testoutput = test_ts
    else #consider training data
        testinput = train_us
        testoutput = train_ts
    end
    ts = [coords[2] for coords in argmax(testoutput, dims=2)] #hardmax function
    ntest, nus = size(testinput)
    correct = 0
    cost = 0
    results = zeros(10, 10) #store prediction results
    for i = 1:ntest #for each datapoint
        target = ts[i]
        us = reshape(testinput[i,:], 1, nus) #input data
        if bias us = [us 1] end #add bias unit
        vs, zs = calc_vs_zs(us, w1s=w1s, w2s=w2s, bias=bias, transfer=transfer)
        if loss #calculate loss function
            cost += sum((zs.*(1.08/sum(zs)) .- testoutput[i,:]).^2)
        end
        prediction = argmax(zs[:]) #use max function for prediction
        results[target, prediction] += 1
        if prediction == target correct += 1 end #count number of correct predictions
    end
    if loss return log10(cost) end
    frac_corr = correct / ntest #proportion of correct predictions

    if Plot
        #plot prediction results
        results ./= sum(results, dims=2)
        println(results)
        println("performance: ", frac_corr)
        #plot heatmap of cross-prediction
        heatmap(results, 0:9, 0:9; xlab="predicted", ylab="target",
            filename="figures/"*filename*"_predictions.png", Title="",
            vmin=0, vmax=0.05)

        figure(figsize=(5,4)) #plot barplot of performance
        bar(0:9, [results[i,i] for i =1:10])
        xticks(0:9, [string(lab) for lab in 0:9])
        xlabel("target")
        ylabel("performance")
        hlines([frac_corr],-2,12, "k","--")
        xlim(-1,10)
        savefig("figures/"*filename*"_predictionhist.png", dpi=120)
        close()
    end
    return frac_corr #return prediction performance
end

function train(; nepoch = 100, fname = "test", w1s=w1s, w2s=w2s, eta=0.1, N=N, Plot=true,
                bias=false, return_weights = false, dropout=[0 0], transfer="tanh",
                train_ts=train_ts, train_us=train_us)
    #trains the network for nepoch epochs and plots a graph of performance'''
    println("training a network with "*string(N)*
            " hidden units and learning rate "*string(eta)*
            " for "*string(nepoch)*" epochs. bias "*string(bias)*
            ", transfer "*transfer)
    n = 0 #number of epochs
    performance = zeros(nepoch+1) #store test performance
    frac_corr0 = test_all(w1s=w1s, w2s=w2s, bias=bias, transfer=transfer)
    println("now at n="*string(n)*
            " performance="*string(round(frac_corr0, digits=2)))
    performance[1] = frac_corr0

    performance_train = zeros(nepoch+1) #also calculate training performance
    performance_train[1] = test_all(w1s=w1s, w2s=w2s, bias=bias, test=false, transfer=transfer)

    while n < nepoch #for each epoch
        n += 1
        oldw1s = copy(w1s)
        oldw2s = copy(w2s)
        t = time() #also calculate time per epoch
        w1s, w2s = train_round(w1s=w1s, w2s=w2s, eta=eta, N=N, bias=bias, dropout=dropout,
                                transfer=transfer) #train a single epoch
        t = time() - t
        err = norm(w1s-oldw1s)+norm(w2s-oldw2s) #how much did weights change?
        #quantify how good we are
        frac_corr = test_all(w1s=w1s, w2s=w2s, bias=bias, loss=false, transfer=transfer)
        frac_train = test_all(w1s=w1s, w2s=w2s, bias=bias, test=false, loss=false, transfer=transfer)
        println("now at n="*string(n)*" test="*string(round(frac_corr, digits=3))*
                " train="*string(round(frac_train, digits=3))*
                " time="*string(round(t, digits=2)))
        performance[n+1] = frac_corr #store results
        performance_train[n+1] = frac_train
    end

    if Plot
        #plot test and training performance over epochs
        mypickle(performance, "pickled_data/"*fname*"_performance.pickled")
        figure(figsize = (5,3))
        plot(0:nepoch, performance)
        plot(0:nepoch, performance_train, "0.5")
        xlabel("epoch")
        ylabel("performance")
        xlim(0, nepoch)
        if size(train_us)[2] == 10 ylim(0.45, 1)
        else ylim(0,1) end
        legend(["test", "train"])
        savefig("figures/"*fname*"_performance.png", dpi=120, bbox_inches="tight")
        close()
    end
    if return_weights return performance, w1s, w2s end
    return performance
end

function repeat_performance(;N = N, ntrial = 10, nepoch = 100, fname = "test",
                    w1s=w1s, w2s=w2s, eta=0.1, Plot = true, bias=false,
                    train_us=train_us, transfer="tanh", dropout = [0 0])
    #train ntrial nets and quantify performance
    performances = zeros(ntrial, nepoch+1)
    for i = 1:ntrial #for each trial
        println("\n\nNew trial: "*string(i))
        ninput, nus = size(train_us)
        nzs = size(train_ts)[2]
        d = Normal(0, 0.01)
        if bias
            w1s = rand(d, nus+1, N) #random initial weights
            w2s = rand(d, N+1, nzs) #random initial weights
        else
            w1s = rand(d, nus, N) #random initial weights
            w2s = rand(d, N, nzs) #random initial weights
        end
        #find learning curve
        performance = train(; eta=eta, nepoch = nepoch, fname = "test",
                            w1s=w1s, w2s=w2s, N=N, Plot=false, bias=bias,
                            transfer=transfer, dropout = dropout)
        performances[i, :] = performance #store performance
    end
    performance = mean(performances, dims = 1)[:] #report mean performance
    return performance
end

function repeat_train(;N = N, ntrial = 10, nepoch = 100, fname = "test",
                    w1s=w1s, w2s=w2s, eta=0.1, Plot = true, bias=false,
                    train_us = train_us, train_ts = train_ts, transfer="tanh",
                    dropout = [0 0])
    #repeat training and return only maximum performance
    performances = zeros(ntrial, nepoch+1)
    for i = 1:ntrial #for each trial
        println("\n\nNew trial: "*string(i))
        ninput, nus = size(train_us)
        nzs = size(train_ts)[2]
        d = Normal(0, 0.01)
        if bias
            w1s = rand(d, nus+1, N) #random initial weights
            w2s = rand(d, N+1, nzs) #random initial weights
        else
            w1s = rand(d, nus, N) #random initial weights
            w2s = rand(d, N, nzs) #random initial weights
        end
        #calculate learning curve
        performance = train(; eta=eta, nepoch = nepoch, fname = "test",
                            w1s=w1s, w2s=w2s, N=N, Plot=false, bias=bias,
                            transfer=transfer, dropout = dropout)
        performances[i, :] = performance
    end
    performance = mean(performances, dims = 1)[:] #find mean

    if Plot #plot result
        figure(figsize = (4,2))
        #figure(figsize = (5,3))
        plot(0:nepoch, performance)
        xlabel("epoch")
        ylabel("performance")
        xlim(0, nepoch)
        if size(train_us)[2] == 10 ylim(0.45, 1)
        else ylim(0,1) end
        savefig("figures/"*fname*"_performance.png", dpi=120, bbox_inches="tight")
        close()
    end

    maxperf, maxepoch = findmax(performance) #find maximum performance
    maxepoch -= 1 #discount zero epoch performance
    return maxperf, maxepoch #return maximum performance and when it occurred
end

function plot_image(us; fname="test")
    #plot our input as a greyscale image
    image = Array(transpose(reshape(us, 28, 28)))
    heatmap(image, [], []; xlab="", ylab="",
    filename="figures/"*fname*".png", Title="")
end

function plot_weights(ws; fname="test")
    #plot weights to our hidden layer
    n1, n2 = 2, 5
    fig, axes = subplots(n1, n2, figsize=(n2,n1))

    for i = 1:(n1*n2)#size(ws)[2]
        println(i)
        image = Array(transpose(reshape(ws[1:784,i], 28, 28)))
        ax = axes[Int(floor((i-1)/n2)+1), (i-1)%n2+1 ]
        ax[:imshow](image, cmap=ColorMap("gray_r"), vmin = minimum(image),
                vmax = maximum(image))
        ax[:set_xticks]([], [])
        ax[:set_yticks]([], [])
    end
    savefig("figures/"*fname*"_weightfig.png", dpi=480, bbox_inches="tight")
    close()
end

data = unpickle("MNIST/data.pickled") #unpack our training and test data
train_us, train_ts, test_us, test_ts = data[1], data[5], data[2], data[6]

parity = false #if parity, work on the parity task!
if parity
    n1 = 400
    n2 = 1000
    train_us = rand([1, -1], (n1, 10)) #random training data
    train_ts =  zeros(n1, 2) #find real answers
    ts = prod(train_us, dims=2)
    train_ts[ts[:,1].==1, 1] .= 1
    train_ts[ts[:,1].==-1, 2] .= 1

    test_us = rand([1,-1], (n2, 10)) #random test data
    test_ts =  zeros(n2, 2) #find real data
    ts = prod(test_us, dims=2)
    test_ts[ts[:,1].==1, 1] .= 1
    test_ts[ts[:,1].==-1, 2] .= 1
end

subsample = false #subsample our training data for overfitting
if subsample
    inds = rand(50000:60000, 1000)
    test_us = train_us[inds, :] #random test data
    test_ts = train_ts[inds, :]
    inds = rand(1:50000, 200)
    train_us = train_us[inds, :] #random training data
    train_ts = train_ts[inds, :]
end

shift = false #if shift, shift out data to zero mean
if shift
    train_us = train_us .- mean(train_us, dims=2)
    test_us = test_us .- mean(test_us, dims=2)
end

ninput, nus = size(train_us) #number of inputs and input units
ntest = size(test_us)[1] #store number of test datapoints
nzs = size(train_ts)[2]

N = 10 #number of hidden units
eta = 0.01 #learning rate
nepoch = 20 #number of epochs to train for
bias = false #do we include bias units?
transfer = "tanh" #specify transfer function (tanh or sigmoid)

d = Normal(0, 0.01) #start with random weights in linear regime
if bias
    w1s = rand(d, nus+1, N) #random initial weights
    w2s = rand(d, N+1, nzs) #random initial weights
else
    w1s = rand(d, nus, N) #random initial weights
    w2s = rand(d, N, nzs) #random initial weights
end

fname="test_N"*string(N)*
        "_eta"*string(eta)[3:end]*
        "_nepoch"*string(nepoch)*
        "_transfer_"*transfer
if parity fname = fname*"_parity" end
#fname = fname*"_nobias"

#performance, w1s, w2s = train(eta=eta, nepoch=nepoch, fname=fname, bias=bias,
#                        return_weights = true, transfer=transfer)

#writedlm("w1s_test_tanh.dlm", w1s) #store out weights!
#writedlm("w2s_test_tanh.dlm", w2s)



ntrial = 10
fname = "test_N"*string(N)*
        "_eta"*string(eta)[3:end]*
        "_nepoch"*string(nepoch)*
        "_ntrial"*string(ntrial)*
        "_transfer_"*transfer
if parity fname = fname*"_parity" end
#repeat_train(eta=eta, nepoch=nepoch, N=N, ntrial=ntrial,
#            fname=fname, bias=bias, dropout=[0.0 0.0], transfer="tanh")


#w1s = readdlm("w1s_N10_tanh.dlm")
#w2s = readdlm("w2s_N10_tanh.dlm")
#test_all(bias=true, Plot=true)
#plot_weights(w1s, fname="w1s_N10")


 #w3s = w1s[1:784,:]*w2s[1:size(w1s)[2],:]
 #plot_weights(w3s, fname="w3s_N300")
