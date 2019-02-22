using Random, Distributions, PyCall, PyPlot, LinearAlgebra, DelimitedFiles

state = []

function move_to_grid(state0, xlo, xhi, ylo, yhi)
    newstate = copy(state0)
    r, v = newstate
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
    newstate[1] = r1
    return newstate
end

function advance(state0, t0, t1)
    #compute the state after time t ignoring collisions
    #update velocity and position
    #velocity v is constant; r1 = v*(t1-t0) + r
    #TESTED
    newstate = copy(state0)
    r, v = newstate
    r1 = r .+ v*(t1-t0)

    ##move back onto grid if we have left grid
    #while r1[1] < xlo
    #    r1[1] += (xhi-xlo)
    #end
    #while r1[1] > xhi
    #    r1[1] -= (xhi-xlo)
    #end
    #while r1[2] < ylo
    #    r1[2] += (yhi-ylo)
    #end
    #while r1[2] > yhi
    #    r1[2] -= (yhi-ylo)
    #end

    newstate[1] = r1
    return newstate
end

function interaction_time1(state, time, k)
    #take a current time 'time' and a state,
    #return the nearest time of crossing boundary k
    #k = 1,2,3,4 ccw from right boundary
    #TESTED
    r, v = state
    if k == 1 || k == 3
        ind = 1 #x coordinates
    else
        ind = 2 #y coordinates
    end
    if k == 1 || k ==2
        lim = 400 #right or top edge
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

function interaction_time2(state1, time1, state2, time2, a0; Print = false)
    #TESTED

    tstar = max(time1, time2)
    v = state2[2] - state1[2]
    r20 = state2[1] + state2[2] * (tstar-time2)
    Print && println("v: ", v)
    sects = [[0,0],[0,400],[0,-400],[400,0],[400,400],[400,-400],[-400,0],[-400,400],[-400,-400]]
    is = []
    times = []
    for i in 1:9
        r10 = state1[1] + state1[2] * (tstar-time1) + sects[i]

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
                if deltat > -0.001 #need event to happen in the future but with error
                    Print && println("new interaction, mods: ", sects[i], " time:", tstar+deltat)
                    times = vcat(times, tstar+max(0, deltat) )
                    is = vcat(is, i)
                end
            end
        end
    end

    Print && println("times: ", times)


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
        newval = 400 #jump to right or top
    end
    Print && println(ind ," ", newval)
    newstate = copy(state)
    newstate[1][ind] = newval
    return newstate
end

function jump2(state1, state2, a0, sectorind, Print = false)
    #process collision between sphere 1 and 2
    #TESTED
    #need to take into account 8 adjacent sectors as well
    sects = [[0,0],[0,400],[0,-400],[400,0],[400,400],[400,-400],[-400,0],[-400,400],[-400,-400]]

    h = a0 #h = a'(t); a(t) = a0*t

    r1, v1 = copy(state1); r2, v2 = copy(state2)
    r1 = r1 + sects[sectorind] #put to correct sector
    newstate1 = copy(state1); newstate2 = copy(state2)

    r12 = r1-r2; r12 = r12/norm(r12)
    Print && println("r12: ", r12)
    v1p = dot(v1, r12) * r12; v1t = v1-v1p
    v2p = dot(v2, r12) * r12; v2t = v2-v2p
    Print && println("v1p: ", v1p, "  v1t: ", v1t)

    v1s = (v2p + h*r12) + v1t #transverse velocity unchanged
    v2s = (v1p - h*r12) + v2t #parallel swapped w/ additive h*(+/-)u12

    newstate1[2] = v1s
    newstate2[2] = v2s

    return newstate1, newstate2
end

function initialize(n, xlo, xhi, ylo, yhi, Print = false)
    #initialize events
    partner = zeros(n); time = zeros(n) #row for disk, column for old vs new
    state = reshape(zeros(2), 1, 2)
    for i = 1:n
        partner[i, 1] = NaN #no new partners
        r = rand(2,1).*[xhi-xlo, yhi-ylo]+[xlo, ylo] #random positioning
        v = ( rand(2,1) .- 0.5 )*2 #initial velocities between -1 and 1
        Print && println( hcat([[r, v]], [[r, v]]), "\n", state )
        state = vcat( state, hcat([ [ r, v  ] ], [ [ r, v ] ] ) )
    end
    state = state[2:end, :]
    time = hcat(copy(time), copy(time))
    partner = hcat(copy(partner), copy(partner))
    #state = hcat(copy(state), copy(state))
    Print && println("copied?")
    return time, partner, state
end

function get_interactions(i)

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
    indi_prev, indj_prev = 0, 0

    while (current_time <= end_time) & (niter < maxiter)
        niter += 1
        Print && println("\n\nNew iteration ", niter)
        if (niter % 10 == 0)
            println("new n: ", niter, "  current size: ", current_time*a0)
            if niter % 500 == 0
                #for i in 1:n
                #    println(state[i, Old[i]][2])
                #    state[i,Old[i]][2] = [0.0; 0.0]
                #end
            end
        end
        #current_time, indi = findmin(time[:, New[i]) #nearest new event; ind is index
        #current_time = minimum( vcat(time[:,1][ New .== 1 ], time[:,2][ New .== 2 ]))
        newtimes = zeros(n)
        #Print && println(New)
        newtimes[New .== 1] = time[:,1][ New .== 1 ]
        newtimes[New .== 2] = time[:,2][ New .== 2 ]
        current_time, indi = findmin(newtimes)
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
        for j in indi_bar
            tcol, sectorind = interaction_time2(
                state[indi, Int(Old[indi])], time[indi, Int(Old[indi])],
                state[j, Int(Old[j])], time[j, Int(Old[j])], a0)
            Ps[j] = tcol
            Print && println("j: ",j,"  PrevPartnerJ: ", partner[j, Old[j]], "  PrevPartneri: ", partner[indi, Old[indi]],
                            "  CurrPartnerJ: ", partner[j, New[j]], "  CurrPartneri: ", partner[indi, New[indi]])
            if partner[j, Old[j]] == indi && partner[indi, Old[indi]] == j && abs(time[indi, Old[indi]]-tcol) < 0.001
                Ps[j] = Inf #already next scheduled event
            end
            sectorinds[j] = sectorind
            #timebar[j] = time[j, New[j]]
        end
        if indi == indi_prev #don't process same collision twice in a row
            Ps[indj_prev] = Inf
        #elseif indi == indj_prev
        #    Ps[indi_prev] = Inf
        end

        Print && println("Ps: ", Ps, "  indi_bar: ", indi_bar)
        Print && println("Sectorinds ", sectorinds)
        Print && println("newtimes: ", newtimes )
        A_indi = (1:n)[ newtimes .>= Ps ]
        Print && println("A_indi: ", A_indi)
        #P is next interaction time of two spheres
        if length(A_indi) > 0
            Pmin, indj = findmin(Ps[A_indi])  #j = A(ind), A(ind) = {j | 1 <= j <= N, j !=ind, time[j, new[j]] >= P[]ind, j]}
            indj = A_indi[indj] #get j value corresponding to this minimum
            sectorind = Int(sectorinds[indj])
            Print && println("Pmin: ", Pmin, "  indj: ", indj)
            indi_prev = indi
            indj_prev = indj
        else
            Pmin, indj = Inf, 0
            indi_prev = 0 #process boundary crossing, reset these
            indj_prev = 0
        end

        Qs = zeros(4)
        for k in 1:4
            Qs[k] = interaction_time1(
                state[indi, Int(Old[indi])], time[indi, Int(Old[indi])], k)
        end
        Print && println("Qs: ", Qs)
        Qmin, indk = findmin(Qs)
        #Print && println("Qmin: ", Qmin, "  indk: ", indk)

        R = min(Pmin, Qmin)
        #Print && println("R: ", R)
        time[indi, Int(New[indi])] = R

        if R < Inf
            #Print && println("time0: ",  time[indi, Old[indi]])
            state1 = advance(state[indi, Old[indi]], time[indi, Old[indi]], R)
            Print && println("advanced")
            Print && println("state1: ", state1)
            if Qmin < Pmin
                state[indi, Int(New[indi])] = jump1(state1, indk)
                Print && println("jumped: ", state[indi, Int(New[indi])])
                partner[indi, Int(New[indi])] = NaN
            else
                time[indj, Int(New[indj])] = R
                state2 = advance(state[indj, Old[indj]], time[indj, Old[indj]], R)
                Print && println("advanced")
                Print && println("state2: ", state2)
                dist = norm(state2[1]-(state1[1]+sects[sectorind]))
                if round(dist) != round(a0*R)
                    println("ERROR! dists do not match. Dist: ", dist, "  a(t): ", a0*R)
                end
                Print && println("Distance: ", norm(state2[1]-(state1[1]+sects[sectorind])), "  a=", a0*R)
                new1, new2 = jump2(state1, state2, a0, sectorind)
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

function max_growing_disks(;thresh=10)
    #thresh is min disk ratio; see how many points we can add
    #while still having radius > thresh
end

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
    PyPlot.savefig(filename)
    close()

end

#function plot_disks(states, endsize)
#    N = size(states)[1]
#    points = zeros(N,2)
#    for i in 1:N
#        points[i, :] = states[i][1]
#    end
#    println(points)
#    plot_points(points, pointsize = endsize/2)
#    return points
#end



function interaction_time2_old(state1, time1, state2, time2; a0=0.2, Print = false)
    #TESTED

    tstar = max(time1, time2)
    v = state2[2] - state1[2]
    r20 = state2[1] + state2[2] * (tstar-time2)
    Print && println("v: ", v)
    sects = [[0,0],[0,400],[0,-400],[400,0],[400,400],[400,-400],[-400,0],[-400,400],[-400,-400]]
    is = []
    times = []
    for i in 1:9
        r10 = state1[1] + state1[2] * (tstar-time1) + sects[i]

        Print && println("r20: ", r20)
        Print && println("r10: ", r10)
        r = r20-r10
        Print && println("r: ", r)

        A = norm(v)^2 - a0^2
        Print && println("A: ", A)
        B = dot(r, v) - a0 * a0*tstar
        Print && println("B: ", B)
        C = norm(r)^2 - (a0*tstar)^2
        Print && println("C: ", C)

        if (B^2 - A*C) > 0
            t = ( (-B - sqrt(B^2 - A*C)) / A )
            if t > tstar #need event to happen in the future
                times = vcat(times, t)
                is = vcat(is, i)
            end
        end
    end

    println("times: ", times)


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
