#functions for plotting grids, heatmaps and other graphs
function plot_points(points; pointsize=10, xlo=0, xhi=400, ylo=0, yhi=400,
            xlab="", ylab="", title="", filename="default.png")
    #given a list of points and some plotting parameters,
    #plots the points. Use PyPlot backend. pointsize is the radius of points

    figsize = [5,5] .* [yhi-ylo, xhi-xlo] ./ mean([yhi-ylo, xhi-xlo])
    fig, ax = subplots(figsize = figsize)
    xlabel(xlab)
    ylabel(ylab)
    xlim([xlo,xhi])
    ylim([ylo,yhi])
    circles = []
    for i in 1:size(points)[1]
        #plot points as circles with radius
        circles = vcat(circles,
        matplotlib[:patches][:Circle](
        (points[i,1], points[i,2]), radius = pointsize, facecolor="b"))
    end
    p = matplotlib[:collections][:PatchCollection](circles)
    ax[:add_collection](p) #plot all points
    println("plotted")
    savefig(filename)
    close()
end

function heatmap(results, xs, ys; xlab="", ylab="",
    filename="default", Title="default")
    #given a matrix of z values and lists of x,y values
    #plots a heatmap

    figure()
    imshow(results, cmap=ColorMap("gray_r"), vmin = minimum(results),
            vmax = maximum(results))
    print("plotted")
    xticks(0:(length(xs)-1), xs, rotation=90)
    yticks(0:(length(ys)-1), ys)
    colorbar()
    xlabel(xlab)
    ylabel(ylab)
    title(Title)
    savefig(filename, bbox_inches = "tight")
    close()
end

function plot_descent()
    #reads result of parameter space investigation and steepest descent_results
    #and projects the descent onto a 3d plot of similarity values

    scan = readdlm("scan_results.txt")
    N = Int(sqrt(size(scan)[1]))
    ms = reshape(scan[:,1], N, N)
    ss = reshape(scan[:,2], N, N)
    us = reshape(scan[:,3], N, N)

    #first plot heatmap of parameter space
    heatmap(convert(Array, us'), 0:0.3:6.3, 0:1:21, ylab="mean", xlab="standard deviation",
            filename="similarity_gray.png", Title="variation of u-score")

    descent = readdlm("descent_results.txt")
    figure()
    #surface plot of parameter space
    plot_surface(ms, ss, us, cmap=ColorMap("gray_r") )
    xlabel("mean")
    ylabel("std")
    zlabel("u1")
    # add 0.03 to descent to raise it above surface such that it's visible
    p = plot(descent[:,1], descent[:,2],
        descent[:,3].+0.3, "k-", linewidth=2)


    savefig("projected_descent.png", bbox_inches="tight")


    close()

    figure(figsize = (6,4))
    #surface plot of parameter space
    plot(descent[:,3], "b-" )
    # add 0.03 to descent to raise it above surface such that it's visible
    xlabel("interation number")
    ylabel("u1")

    savefig("descent.png", bbox_inches="tight")
    close()
end


function plot_packings()
    #plot square and hexagonal packings on 400x400 grid for r=10 spheres
    xs=repeat(0:20:400, outer=21) #square packing x values
    ys=repeat(0:20:400, inner=21) #square packing y values
    plot_points(hcat(xs, ys), filename = "square.png", pointsize=10,
                xlo = 0, xhi = 400, ylo = 0, yhi = 400)

    #define hexagonal x and y parameters
    xs = vcat(repeat(0:20:400, outer=12), repeat(10:20:390, outer=12))
    ys = vcat(repeat(0:2*sqrt(300):400, inner=21),
        repeat(sqrt(300):2*sqrt(300):400, inner=20))
    plot_points(hcat(xs, ys), filename = "hexagonal.png", pointsize=10,
                xlo = 0, xhi = 400, ylo = 0, yhi = 400)
end
