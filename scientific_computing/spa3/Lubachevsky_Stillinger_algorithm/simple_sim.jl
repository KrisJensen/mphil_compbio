using Random, Distributions, PyCall, PyPlot, LinearAlgebra, DelimitedFiles

state = []
xhi = 420; yhi = 420 #specify default gridsize with global variables

function move_to_grid(state0, xlo, xhi, ylo, yhi)
    #given a point in a state, moves it onto the grid by
    #applying periodic boundary conditions
    newstate = copy(state0)
    r = newstate[1:2]; v = newstate[3:4]
    r1 = copy(r)

    #do this 'naively' by considering all directions separately
    while r1[1] < xlo
        r1[1] += (xhi-xlo)
    end
    while r1[1] > xhi
        r1[1] -= (xhi-xlo)
    end
    while r1[2] < ylo
        r1[2] += (yhi-ylo)
    end
    while r1[2] > yhi
        r1[2] -= (yhi-ylo)
    end
    newstate[1:2] = r1
    return newstate
end

function advance(state0, dt)
    #compute the state after time t ignoring collisions
    #velocity v is constant; r1 = v*(t1-t0) + r
    newstate = copy(state0)
    r=newstate[1:2]; v = newstate[3:4]
    r1 = r .+ v*dt #find new position
    newstate[1:2] = r1
    return newstate
end

function interaction_time1(state, time, k)
    #take a current time 'time' and a state,
    #return the nearest time of crossing boundary k
    #k = 1,2,3,4 ccw from right boundary
    r = state[1:2]; v = state[3:4]
    if k == 1 || k == 3
        ind = 1 #x coordinates
    else
        ind = 2 #y coordinates
    end
    if k == 1 || k ==2
        lim = xhi #right or top edge
    else
        lim = 0 #bottom edge
    end
    t = 1/v[ind] * ( lim + v[ind] * time - r[ind] ) #time to crossing
    if t > time
        return t #if happens in future that's fine
    else
        return Inf #not moving towards boundary, return Inf
    end

end

function interaction_time2(state1, state2, tstar, a0; Print = false)
    #given two states, return a list of times at which they collide
    #consider main grid and one period in either direction

    v = state2[3:4] - state1[3:4] #difference in velocity
    r20 = state2[1:2]
    Print && println("v: ", v)
    sects = [[0,0],[0,xhi],[0,-xhi],[xhi,0],[xhi,xhi],[xhi,-xhi],[-xhi,0],[-xhi,xhi],[-xhi,-xhi]]
    is = [] #store sector indices for collisions
    times = []
    for i in 1:9
        r10 = state1[1:2] + sects[i] #apply period boundary conditions to state1
        Print && println("r20: ", r20)
        Print && println("r10: ", r10)
        r = r20-r10 #difference in position
        Print && println("r: ", r)

        #solve |r+vt|^2=[a(t*)+aot]^2
        A = norm(v)^2 - a0^2
        Print && println("A: ", A)
        B = 2 * dot(r, v) - 2*a0^2*tstar
        Print && println("B: ", B)
        C = norm(r)^2-a0^2*tstar^2
        Print && println("C: ", C)

        if (B^2 - 4*A*C) > 0 #quadratic has solutions
            deltat1 = ( (-B + sqrt(B^2 - 4*A*C)) / (2*A) )
            deltat2 = ( (-B - sqrt(B^2 - 4*A*C)) / (2*A) )
            for deltat in [deltat1, deltat2]
                if deltat > -0.0001 #need event to happen in the futureallowing for numerical errors
                    Print && println("new interaction, mods: ", sects[i], " time:", tstar+deltat)
                    times = vcat(times, tstar+max(0, deltat) ) #append time of collision
                    is = vcat(is, i) #note which sector state1 is in
                end
            end
        end
    end
    Print && println("times: ", times)
    return times, is #return list of all collisions
end

function jump1(state, k, Print = false)
    #jump across edge. Given state and index of edge, returns new state
    if k == 1 || k == 3
        ind = 1 #x coordinates
    else
        ind = 2 #y coordinates
    end
    if k == 1 || k ==2
        newval = 0 #right or top edge; jump to left or bottom
    else
        newval = xhi #jump to right or top
    end
    Print && println(ind ," ", newval)
    newstate = copy(state)
    newstate[ind] = newval #update position
    return newstate
end

function jump2(state1, state2, a0, sectorind; Print = false)
    #process collision between sphere 1 and 2. Update velocities
    #need to provide information on which sector state2 is in
    sects = [[0,0],[0,xhi],[0,-xhi],[xhi,0],[xhi,xhi],[xhi,-xhi],[-xhi,0],[-xhi,xhi],[-xhi,-xhi]]
    h = a0 #h = a'(t); a(t) = a0*t

    r1, v1 = copy(state1[1:2]), copy(state1[3:4])
    r2, v2 = copy(state2[1:2]), copy(state2[3:4])
    r1 = r1 + sects[sectorind] #place in correct sector
    newstate1 = copy(state1); newstate2 = copy(state2)

    r12 = r1-r2; r12 = r12/norm(r12) #normalized inter-point vector
    Print && println("r12: ", r12)
    v1p = dot(v1, r12) * r12; v1t = v1-v1p #compose v into perpendicular and transverse components
    v2p = dot(v2, r12) * r12; v2t = v2-v2p
    Print && println("v1p: ", v1p, "  v1t: ", v1t)

    v1s = (v2p + h*r12) + v1t #transverse velocity unchanged
    v2s = (v1p - h*r12) + v2t #parallel swapped w/ additive h*(+/-)u12

    newstate1[3:4] = v1s #update velocities
    newstate2[3:4] = v2s
    return newstate1, newstate2
end

function initialize(n, xlo, xhi, ylo, yhi, Print = false)
    #initialize array of positions and velocities
    state = zeros(n, 4) #store x,y position and velocity
    for i = 1:n
        r = rand(2,1).*[xhi-xlo, yhi-ylo]+[xlo, ylo] #random positioning
        v = ( rand(2,1) .- 0.5 )*2 #initial velocities between -1 and 1 in each direction
        Print && println( hcat([[r, v]], [[r, v]]), "\n", state )
        state[i, 1:2] = r; state[i, 3:4] = v
    end
    return state
end

function LS_simulation(;n=10, a0 = 1.00, xlo=0, xhi=xhi, ylo=0, yhi=yhi, end_time=1000000, maxiter = 100000,
                    Print = false, Printall = false, nhist = 100, writepoints = 20000, bisection = false)
    #naive implementation scales as n^2 but given the relatively small system sizes we're working with,
    #this runs in a reasonable amount of time and is much simpler to implement. Hence
    #the preferred implementation for the present purposes.
    #a(t) = a0*t #DIAMETER GROWTH RATE
    #h = a'(t) = a0

    nroll = 0
    sects = [[0,0],[0,xhi],[0,-xhi],[xhi,0],[xhi,xhi],[xhi,-xhi],[-xhi,0],[-xhi,xhi],[-xhi,-xhi]]
    niter = 0 #number of iterations run
    current_time = 0 #start at time zero
    state = initialize(n, xlo, xhi, ylo, xhi)
    recent_events = zeros(nhist, 3) #store nhist most recent events as [i, partner, time]
    recount = 0 #keep track of where in the recent events array we are

    while (current_time <= end_time) & (niter < maxiter)
        oldstate = copy(state) #keep a copy of the old state in case we need to reverse iteration
        recount = (recount % nhist) + 1 #increment by 1
        niter += 1
        mvel = mean(abs.(state[:,3:4]))
        if mvel > 10 #reduce velocities if too high; prevent divergence
            state[:,3:4] /= (4*mvel)
            recent_events = zeros(nhist, 3)
            println("resat velocities")
        end
        Print && println("\n\nNew iteration ", niter)
        if (niter % 10 == 0)
            println("new n: ", niter, "  current size: ", current_time*a0)
            #print progress
            if niter % writepoints == 0
                #occasionally store positions of all points
                writedlm("intermediate_points/points_"*string(n)*"_"*string(niter)*
                            "_"*string(round(current_time*a0,digits=2)), state[:,1:2])
            end
        end
        Print && for i in 1:n println(i, ": ", state[i,:]) end
        collision = true; tmin = Inf; imin=NaN; partner=NaN; kmin=NaN; sectind = NaN #reset parameters
        for i in 1:n
            for j in (i+1):n
                #consider all possible pairwise collisions; find interaction times
                tcols, sectorinds = interaction_time2(state[i,:], state[j,:], current_time, a0)
                for ind in 1:length(tcols)
                    #consider collision from every sector
                    tcol = tcols[ind]; sectorind = sectorinds[ind]
                    if tcol < tmin
                        #new next event
                        add = true
                        for prev in 1:nhist
                            if i == recent_events[prev,1] && j == recent_events[prev,2] && abs(tcol - recent_events[prev,3]) < 0.0001
                                add = false #already processed this event
                                Print && println("might have processed this event? ", recent_events[prev,:])
                            end
                        end
                        if add
                            #store parameters
                            collision = true
                            imin = i
                            tmin = tcol
                            partner = j
                            sectind = sectorind
                            Print && println("new tmin ", tmin, " collide ", imin, " with ", partner, " sector: ", sectind)
                        end
                    end
                end
            end

            for k in 1:4
                #consider collisions with every edge
                tcross = interaction_time1(state[i,:], current_time, k)
                if tcross < tmin
                    add = true
                    for prev in 1:nhist
                        if i == recent_events[prev,1] && isnan(recent_events[prev,2]) && abs(tcross - recent_events[prev,3]) < 0.000001
                            add = false #already processed this event
                            Print && println("might have processed this event? ", recent_events[prev,:])
                        end
                    end
                    if add
                        #if new next event, store parameters
                        imin = i
                        collision = false
                        tmin = tcross
                        kmin = k
                        Print && println("new tmin ", tmin, " cross ", kmin,)
                    end

                end
            end
        end

        for i in 1:n
            state[i,:] = advance(state[i,:], tmin-current_time)
            #let simulation progress till next event
        end

        if collision
            Print && println("colliding ", imin, " with ", partner, " at ", tmin)
            recent_events[recount, :] = [imin, partner, tmin] #store event
            dist = norm(state[partner,1:2]-(state[imin, 1:2]+sects[sectind]))
            if round(dist) != round(a0*tmin)
                #check that the disks are in fact colliding
                println("ERROR! dists do not match. Dist: ", dist, "  a(t): ", a0*tmin)
            end
            Print && println("Distance: ", norm(state[partner,1:2]-(state[imin,1:2]+sects[sectind])), "  a=", a0*tmin)
            #finally process collision by updating velocities
            state[imin,:], state[partner,:] = jump2(state[imin,:], state[partner,:], a0, sectind)
        else
            #if next event is a border crossing
            Print && println(imin, " crossing ", kmin, " at ", tmin)
            recent_events[recount, :] = [imin, NaN, tmin]
        end

        for i in 1:n
            #move all our points to the grid using periodic boundary conditions
            state[i,:] = move_to_grid(state[i,:], xlo, xhi, ylo, yhi)
        end

        if !collision
            #jump border after moving to grid, otherwise numerical inaccuracies
            #can lead us to processing the same event over and over again
            state[imin,:] = jump1(state[imin,:], kmin)
            Print && println("jumped: ", state[imin,:])
        end


        dmin = Inf; imin = NaN; jmin = NaN; smin = NaN
        for i in 1:n
            for j in (i+1):n
                dmini = Inf
                for sect in sects
                    #find all nearest neighbor distances
                    d = norm(state[i,1:2] + sect - state[j,1:2])
                    if d < dmini
                        dmini = d
                        if d < dmin
                            dmin = d
                            imin = i
                            jmin = j
                            smin = sect
                        end
                    end
                end
                Printall && println("distance ", i,",", j, ": ", dmini)
            end
        end
        Print && println("dmin: ", dmin, " size: ", tmin*a0, " time: ", tmin)
        #println("dmin: ", dmin, " time: ", tmin)
        if tmin*a0 > (dmin + 0.00001) #if two points are closer than their diameter, something is wrong
            #undo last simulation step and try something new
            println("rolling back event; ", imin, " and ", jmin, " too close in sector ", smin)
            println(state[imin,1:2], state[jmin,1:2])
            state = copy(oldstate) #return state to previous step
            nroll += 1
            if nroll == 50 #if we cannnot do anything without disks overlapping, simulation is saturated
                if bisection return true end
                return state, a0*current_time
            end
        else
            #nothing overlaps; update time and proceed to next event
            nroll = 0
            current_time = tmin
            if bisection
                if tmin*a0 > 20 #for binary search, exit if spheres are larger than 20
                    writedlm("intermediate_points/points_"*string(n)*"_"*string(niter)*
                                "_"*string(round(current_time*a0,digits=2)), state[:,1:2])
                    return false
                end
            end
        end
        current_time = tmin
    end

    if bisection return true end
    return state, a0*current_time
end

#println("\n\n\nnew simulation:")
points, endsize = LS_simulation()
plot_points_LS(points[:,1:2], pointsize = endsize/2)


function plot_points_LS(points; pointsize=10, xlo=0, xhi=xhi, ylo=0, yhi=yhi,
            xlab="", ylab="", title="", filename="default.png")
    #given a list of points and some plotting parameters,
    #plots the points. Plot with periodic boundary condidions

    sects = [[0,0],[0,xhi],[0,-xhi],[xhi,0],[xhi,xhi],[xhi,-xhi],[-xhi,0],[-xhi,xhi],[-xhi,-xhi]]

    figsize = [5,5] .* [yhi-ylo, xhi-xlo] ./ mean([yhi-ylo, xhi-xlo])
    fig, ax = subplots(figsize = figsize)
    xlabel(xlab)
    ylabel(ylab)
    xlim([xlo,xhi])
    ylim([ylo,yhi])
    circles = []
    for i in 1:size(points)[1]
        for sect in sects #plot with periodic boundaries
            #plot points as circles with radius
            circles = vcat(circles,
                        matplotlib[:patches][:Circle](
                        (points[i,1]+sect[1], points[i,2]+sect[2]), radius = pointsize, facecolor="b"))
        end
    end
    p = matplotlib[:collections][:PatchCollection](circles)
    ax[:add_collection](p) #plot all points
    println("plotted")
    PyPlot.savefig(filename)
    close()

end


function bisection_LS()
    nmin, nmax = 1, 520
    #perform bisection search to determine saturation limit for LS model
    while (nmax-nmin) > 1
        global n = Int(round(mean([nmin, nmax]))) #get middle value
        println("new n: ", n, "  nmin: ", nmin, "  nmax: ", nmax, " ")
        global sat = LS_simulation( n = n ) #check if saturated

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

function repeat_max_LS(N = 10)
    #run N bisection searches to find maximum number of points we can add
    #in LS model
    ns = zeros(N)
    for i in 1:N
        ns[i] = bisection_disks() #run bisection search
        println("new n: ", ns[i])
    end
    println("\nn's are:", ns)
    return ns
end

#repeat_max_LS()
