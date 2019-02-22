#functions for determining the maximum number of points
#that can be added to a grid
include("plots.jl")
include("dmin2d.jl")
include("birth_death.jl")

function saturation(n;N=10, method = "seq")
    #try constructing a grid with n points and see if it's saturated
    #consider n to be saturating if we N consecutive times fail to add a
    #single point 10,000 times

    for i in 1:N
        #try N times to add n points
        if method == "seq"
            trial = size(dmin2d(n, 20, 0, 0, 400, 0, 400, Plot=false))[1]
        else
            trial = size(birth_death_model(n, m=20, xlo=0, xhi=400, ylo=0, yhi=400))[1]
        end
        if trial == n
            return false #not saturated yet
        end
    end
    return true #failed N times in a row; saturated
end

function test_max(;nmax = 492, method = "seq")
    #perform bisection search to determine saturation limit
    #do this using the function saturation()
    #initial nmax is theoretical maximum for 400x400 grid of r=10 spheres
    nmin = 0
    while (nmax-nmin) > 1
        global n = Int(round(mean([nmin, nmax]))) #get middle value
        println("new n: ", n, "  nmin: ", nmin, "  nmax: ", nmax, " ")
        global sat = saturation(n, method=method) #check if saturated

        if sat
            nmax = n #new upper limit
            println("saturated")
        else
            nmin = n #new lower limit
            println("unsaturated")
        end
    end
    if sat
        return (n-1) #couldn't add this many points
    else
        return n #could just add this many points
    end
end

function repeat_max(N = 10; method = "seq")
    #run N bisection searches to find maximum number of points we can add
    #this allows us to look at the variation generated by this method
    ns = zeros(N)
    for i in 1:N
        ns[i] = test_max(method = method)
    end
    println("\nn's are:", ns)
    return ns
end