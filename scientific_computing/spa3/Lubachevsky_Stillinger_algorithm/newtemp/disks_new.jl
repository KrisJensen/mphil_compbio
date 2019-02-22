using Random, Distributions, PyCall, PyPlot, LinearAlgebra, DelimitedFiles

state = []

xhi = 410; yhi = 410

function move_to_grid(state0, xlo, xhi, ylo, yhi)
    newstate = copy(state0)
    r = newstate[1:2]; v = newstate[3:4]
    r1 = copy(r)

    #move back onto grid if we have left grid
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
    #update velocity and position
    #velocity v is constant; r1 = v*(t1-t0) + r
    #TESTED
    newstate = copy(state0)
    r=newstate[1:2]; v = newstate[3:4]
    r1 = r .+ v*dt
    newstate[1:2] = r1
    return newstate
end

function interaction_time1(state, time, k)
    #take a current time 'time' and a state,
    #return the nearest time of crossing boundary k
    #k = 1,2,3,4 ccw from right boundary
    #TESTED
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
    #k = 1; xhi = 400
    #400 = r[1]+v[1]*(t-time)
    #400 + v[1]*time = r[1] + v[1]*t
    #t = 1/v[1] * (400 + v[1]*time - r[1])

    t = 1/v[ind] * ( lim + v[ind] * time - r[ind] )

    if t > time
        return t #if happens in future that's fine
    else
        return Inf #happens in past, return Inf
    end

end

function interaction_time2(state1, state2, tstar, a0; Print = false)
    #TESTED

    v = state2[3:4] - state1[3:4]
    r20 = state2[1:2]
    Print && println("v: ", v)
    sects = [[0,0],[0,xhi],[0,-xhi],[xhi,0],[xhi,xhi],[xhi,-xhi],[-xhi,0],[-xhi,xhi],[-xhi,-xhi]]
    is = []
    times = []
    for i in 1:9
        r10 = state1[1:2] + sects[i]

        Print && println("r20: ", r20)
        Print && println("r10: ", r10)
        r = r20-r10
        Print && println("r: ", r)

        A = norm(v)^2 - a0^2
        Print && println("A: ", A)
        B = 2 * dot(r, v) - 2*a0^2*tstar
        Print && println("B: ", B)
        C = norm(r)^2-a0^2*tstar^2
        Print && println("C: ", C)

        if (B^2 - 4*A*C) > 0
            deltat1 = ( (-B + sqrt(B^2 - 4*A*C)) / (2*A) )
            deltat2 = ( (-B - sqrt(B^2 - 4*A*C)) / (2*A) )
            for deltat in [deltat1, deltat2]
                if deltat > -0.0001 #need event to happen in the future but with error
                    Print && println("new interaction, mods: ", sects[i], " time:", tstar+deltat)
                    times = vcat(times, tstar+max(0, deltat) )
                    is = vcat(is, i)
                end
            end
        end
    end

    Print && println("times: ", times)

    return times, is
    if length(times) > 0
        tmin, indmin = findmin(times)
        return tmin, is[indmin]
    else
        return Inf, NaN
    end

    #if ((B > 0) & (A >= 0)) || (B^2-A*C) < 0 || A == 0
    #    return Inf
#    else
#        return (-B-sqrt(B^2-A*C))/A
#    end
end

function jump1(state, k, Print = false)
    #jump across edge. Given state and index of edge, returns new state
    #TESTED
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
    newstate[ind] = newval
    return newstate
end

function jump2(state1, state2, a0, sectorind; Print = false)
    #process collision between sphere 1 and 2
    #TESTED
    #need to take into account 8 adjacent sectors as well
    sects = [[0,0],[0,xhi],[0,-xhi],[xhi,0],[xhi,xhi],[xhi,-xhi],[-xhi,0],[-xhi,xhi],[-xhi,-xhi]]

    h = a0 #h = a'(t); a(t) = a0*t

    r1, v1 = copy(state1[1:2]), copy(state1[3:4])
    r2, v2 = copy(state2[1:2]), copy(state2[3:4])
    r1 = r1 + sects[sectorind] #put to correct sector
    newstate1 = copy(state1); newstate2 = copy(state2)

    r12 = r1-r2; r12 = r12/norm(r12)
    Print && println("r12: ", r12)
    v1p = dot(v1, r12) * r12; v1t = v1-v1p
    v2p = dot(v2, r12) * r12; v2t = v2-v2p
    Print && println("v1p: ", v1p, "  v1t: ", v1t)

    v1s = (v2p + h*r12) + v1t #transverse velocity unchanged
    v2s = (v1p - h*r12) + v2t #parallel swapped w/ additive h*(+/-)u12

    newstate1[3:4] = v1s
    newstate2[3:4] = v2s

    return newstate1, newstate2
end

function initialize(n, xlo, xhi, ylo, yhi, Print = false)
    #initialize events
    state = zeros(1, 2)
    partner = zeros(n); time = zeros(n)
    for i = 1:n
        r = rand(2,1).*[xhi-xlo, yhi-ylo]+[xlo, ylo] #random positioning
        v = ( rand(2,1) .- 0.5 )*2 #initial velocities between -1 and 1
        Print && println( hcat([[r, v]], [[r, v]]), "\n", state )
        state = vcat( state, hcat([ hcat(r, v) ], [ hcat(r, v) ] ) )
        state[i, 1:2] = r; state[i, 3:4] = v
    end
    state = state[2:end, :]
    time = hcat(copy(time), copy(time))
    partner = hcat(copy(partner), copy(partner))
    #state = hcat(copy(state), copy(state))
    Print && println("copied")
    return time, partner, state
end

function growing_disks(;n=100, a0 = 1, xlo=0, xhi=400, ylo=0, yhi=400, end_time=1000000, Print = false, maxiter = 10000)
    #a(t) = a0*t #DIAMETER GROWTH RATE
    #h = a'(t) = a0
    #we store time, state and partner separately rather than making an actual 'event array'
    sects = [[0,0],[0,400],[0,-400],[400,0],[400,400],[400,-400],[-400,0],[-400,400],[-400,-400]]
    niter = 0
    current_time = 0
    New = ones(n); Old = repeat([2], inner=n) #initialize 'current index' for each disk
    time, partner, state = initialize(n, xlo, xhi, ylo, xhi)

    recent_events = zeros(50, 3) #store 10 most recent events as [i, partner, time]
    recount = 0

    while (current_time <= end_time) & (niter < maxiter)
        oldstate = copy(state)
        recount = (recount % 50) + 1
        niter += 1

        Print && println("\n\nNew iteration ", niter)
        if (niter % 10 == 0)
            println("new n: ", niter, "  current size: ", current_time*a0)
        end

        newtimes = zeros(n)
        newtimes[New .== 1] = time[:,1][ New .== 1 ]
        newtimes[New .== 2] = time[:,2][ New .== 2 ]
        current_time, indi = findmin(newtimes)

        vels = 0
        for i in 1:n #update old times
            state[i,Old[i]] = move_to_grid( advance(state[i,Old[i]], current_time - time[i, Old[i]]) )
            vels += mean(abs.(state[i, Old[i]][3:4]))
        end

        mvel = vels/n
        if mvel > 5
            for i in 1:n #update old times
                state[i,Old[i]][3:4] /= (4*mvel)
            end
            recent_events = zeros(50, 3)
            Print && "resat velocities"
        end

        Print && for i in 1:n println(i, ": ", state[i, Old[i]]) end
        #Print && println("newtimes: ", newtimes)
        Print && println("currtime: ", current_time)
        Print && println("indi: ", indi, "  indi_prev: ", indi_prev, "  indj_prev: ", indj_prev)
        New[indi] = Int.(copy(Old[indi]))
        Old[indi] = Int.(3 - copy(New[indi])) #swap new and old
        New = Int.(New); Old = Int.(Old)

        indi_bar = vcat(1:indi-1, indi+1:n) #indices not equal to i
        Ps = zeros(n); Ps[indi] = Inf
        sectorinds = zeros(n)
        #timebar = zeros(n); timebar[indi] = Inf
        collision = true; tmin = Inf; imin=NaN; partner=NaN; kmin=NaN; sectind = NaN
        for j in indi_bar

            tcols, sectorinds = interaction_time2(state[indi, Old[indi]], state[j, Old[j]], current_time, a0)

            for ind in 1:length(tcols)
                tcol = tcols[ind]; sectorind = sectorinds[ind]
                if tcol < tmin && newtimes[j] > tcol
                    add = true
                    for prev in 1:50
                        if i == recent_events[prev,1] && j == recent_events[prev,2] && abs(tcol - recent_events[prev,3]) < 0.0001
                            add = false #already processed this event
                            Print && println("might have processed this event? ", recent_events[prev,:])
                        end
                    end
                    if add
                        collision = true
                        tmin = tcol
                        indj = j
                        sectind = sectorind
                        Print && println("new tmin ", tmin, " collide ", imin, " with ", partner, " sector: ", sectind)
                    end
                end
            end

        end

        for k in 1:4
            tcross = interaction_time1(state[indi,Old[indi]], current_time, k)
            if tcross < tmin
                add = true
                for prev in 1:50
                    if i == recent_events[prev,1] && isnan(recent_events[prev,2]) && abs(tcross - recent_events[prev,3]) < 0.000001
                        add = false #already processed this event
                        Print && println("might have processed this event? ", recent_events[prev,:])
                    end
                end
                if add
                    collision = false
                    tmin = tcross
                    kmin = k
                    Print && println("new tmin ", tmin, " cross ", kmin,)
                end

            end
        end


        time[indi, Int(New[indi])] = tmin

        if R < Inf
            #Print && println("time0: ",  time[indi, Old[indi]])
            state1 = advance(state[indi, Old[indi]], tmin - time[indi, Old[indi]])
            Print && println("advanced")
            Print && println("state1: ", state1)
            if !collision
                state[indi, Int(New[indi])] = jump1(state1, indk)
                Print && println("jumped: ", state[indi, Int(New[indi])])
                partner[indi, Int(New[indi])] = NaN
            else
                time[indj, Int(New[indj])] = tmin
                state2 = advance(state[indj, Old[indj]], tmin - time[indj, Old[indj]])
                Print && println("advanced")
                Print && println("state2: ", state2)
                dist = norm(state2[1]-(state1[1]+sectind))
                if round(dist) != round(a0*tmin)
                    println("ERROR! dists do not match. Dist: ", dist, "  a(t): ", a0*R)
                end
                Print && println("Distance: ", norm(state2[1]-(state1[1]+sectind)), "  a=", a0*tmin)
                new1, new2 = jump2(state1, state2, a0, sectind)
                state[indi, New[indi]] = move_to_grid(new1, xlo, xhi, ylo, xhi)
                state[indj, New[indj]] = move_to_grid(new2, xlo, xhi, ylo, xhi)
                Print && println("jumped: ", state[indi, New[indi]], " ", state[indj, New[indj]])
                indm = partner[indj, New[indj]]
                #Print && println("previous partner2: ", indm)
                partner[indi, New[indi]] = Int(indj)
                partner[indj, New[indj]] = Int(indi)
                if (!isnan(indm)) & (indm != indi) #j has collided with i, is no longer ms next partner
                    indm = Int(indm)
                    state[indm, New[indm]] = advance(state[indm, New[indm]],
                                                    time[indm, Old[indm]],
                                                    time[indm, New[indm]],
                                                    xlo, xhi, ylo, yhi)
                    partner[indm, New[indm]] = NaN
                    #Print && println("advanced ", indm)
                end
            end
        end
    end

    points = zeros(n,2)
    for i in 1:n
        state[i, Old[i]] = move_to_grid(advance(state[i, Old[i]], time[i, Old[i]], current_time), xlo, xhi, ylo, yhi)
        points[i, :] = state[i, Old[i]][1]
    end
    return points, a0*current_time
end

println("\n\n\nnew simulation:")
points, endsize = growing_disks()
plot_points(points, pointsize = endsize/2)

function plot_points(points; pointsize=10, xlo=0, xhi=xhi, ylo=0, yhi=yhi,
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
    PyPlot.savefig(filename)
    close()

end
