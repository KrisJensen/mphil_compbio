include("mlp.jl")

function heatmap(results, xs, ys; xlab="", ylab="",
    filename="default", Title="default")
    #given a matrix of z values and lists of x,y values
    #plots a heatmap
    figure()
    imshow(results, cmap=ColorMap("gray_r"), vmin = minimum(results),
            vmax = maximum(results))
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

function test_params(param1, vals1, param2, vals2; ntrial=5, fname="test", transfer="tanh")
    #define defaults
    nepoch = 15
    N = 30
    eta = 0.1

    n1, n2 = length(vals1), length(vals2)
    performances = zeros(n1, n2)
    epochs = zeros(n1, n2)
    for (i, val1) in enumerate(vals1)
        for (j, val2) in enumerate(vals2)
            if (param1, param2) == ("N", "eta")
                performance, epoch = repeat_train(N=val1, eta=val2, nepoch=nepoch,
                                        ntrial=ntrial, Plot=false, transfer=transfer)
            end
            performances[i,j] = performance
            epochs[i,j] = epoch
            println(param1,": ", val1, "  ",param2,": ",val2,
                    " performance: ", performance, " epoch: ", epoch )
        end
    end

    heatmap(performances, vals2, vals1; xlab=param2, ylab=param1,
        filename="figures/"*fname*"_perf.png", Title=param1*" "*param2) #plot

    heatmap(epochs, vals2, vals1; xlab=param2, ylab=param1,
        filename="figures/"*fname*"_epochs.png", Title=param1*" "*param2) #plot
end

function compare_learning_curves(param, vals;ntrial=5, nepoch=25, fname="test",
                                bias=true, N=50, eta=0.05, transfer="tanh",
                                dropout = [0 0])
    #eta = 0.05
    #N = 40

    figure(figsize = (5,3))
    cols = []
    for (i, val) = enumerate(vals)
        println("new val:", val)
        if param == "N"
            N = val
        elseif param == "eta"
            eta = val
        elseif param == "bias"
            bias = val
        elseif param=="transfer"
            transfer = val
        elseif param=="dropout"
            dropout = val
        end
        performance = repeat_performance(N=N, eta=eta, nepoch=nepoch, dropout=dropout,
                                ntrial=ntrial, Plot=false, bias=bias, transfer=transfer)

        plot(0:nepoch, performance, color=string(1-i/length(vals)))

    end
    xlabel("epoch")
    ylabel("performance")
    legend([string(val) for val in vals])
    savefig("figures/"*fname*".png", dpi=120, bbox_inches="tight")
    close()
end



#test_params("N", [10;20;30;40;50], "eta", [0.01;0.05;0.1;0.3;0.5],
#            fname="test_N_eta")

compare_learning_curves("eta", [0.00001; 0.0001; 0.001; 0.01; 0.1; 1.0], fname="testeta_3vals", N=50, nepoch=15, bias=false, ntrial=3)
#compare_learning_curves("N", [20; 50; 200], fname="testN_3vals")
#compare_learning_curves("bias", [true; false], fname="test_bias_parity",
#                eta=0.075, ntrial=10, nepoch=300, N=50, transfer="tanh")
#compare_learning_curves("transfer", ["tanh"; "sigmoid"], fname="test_transfer", ntrial=5, nepoch=10, N=50)
#compare_learning_curves("dropout", [[0.1 0.1], [0 0]], transfer="sigmoid", fname="test_dropout", ntrial=10, nepoch=300, N=50, eta=0.075)
