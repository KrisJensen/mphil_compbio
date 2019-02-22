#functions for running the 'birth and death' model
include("plots.jl")

function birth_death_round(points, m, xlo, xhi, ylo, yhi)
    #I don't entirely understand what we're meant to do here...
    #pick out points in random order, remove point and try to add subject to constraints
    #hardcode a standard deviation of zero

    n = size(points)[1]
    sequence = randperm(n) #numbers 1:n in random order
    nmax = 10000 #max number points we try to add before concluding it's saturated

    for i in sequence
        ## Point i must now be killed, and a new point
        ## positioned (born) randomly subject to satisfying
        ## the minimal distance constraint.
        niter = 0
        add = false
        while (!add) & (niter < nmax)
        #keep generating new point until it satisfies constraint or we've exceeded our limit
            niter += 1
            global coords = rand(2,1).*[xhi-xlo, yhi-ylo]+[xlo, ylo] #new random coordinates
            add = true
            for j in vcat(1:(i-1), (i+1):n) #check distance to all other points
                if add
                    if norm(coords-points[j,1:2]) < m #too close to point
                        add = false #can't add; try new point
                    end
                end
            end
        end
        if add
            points[i,:] = coords' #satisfies constraint
        else
            println("Have failed to add ", nmax,
                " points in a row. Exiting. Added ", i-1)
            return(zeros(0,0)) #failed to add a point nmax times in a row
        end
    end
    return points
end

function birth_death_model(n; m=20, xlo=0, xhi=400, ylo=0, yhi=400, nepochs = 10)
    #generate a grid of points by doing 10 consecutive rounds of
    #shuffling points subject to dmin. we hardcode here s=0
    d = [m]
    points = zeros(n,2)
    points[:,1] = (rand(n).*(xhi-xlo)).+xlo #random x coordinates
    points[:,2] = (rand(n).*(yhi-ylo)).+ylo #random y coordinates
    for epoch = 1:nepochs
        #shuffle points
        points = birth_death_round(points, m, xlo, xhi, ylo, yhi )
    end
    return points
end
