include("plots.jl")

function find_next_event(state, xlo, xhi, ylo, yhi)

end

function advance(state0, t0, t1, a0)
  #compute the state after time t ignoring collisions
  #update velocity and position
  #velocity v is constant; r1 = v*(t1-t0) + r
  r = state[1], v=state[2]
  r1 = r .+ v*(t1-t0)
  return [r1, v]
end

function interaction_time1(state, time, k, a0)
    #take a current time 'time' and a state,
    #return the nearest time of crossing boundary k
    #k = 1,2,3,4 ccw from right boundary
    r, v = state
    if k == 1 or k == 3
        ind = 1 #x coordinates
    else
        ind = 2 #y coordinates
    end
    if k == 1 or k ==2
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

function interaction_time2(state1, time1, state2, time2; a0=3.2)

    tstar = max(time1, time2)
    v = state2[2] - state1[2]
    r20 = state2[1] + state2[2] * (tstar-time2)
    r10 = state1[1] + state1[2] * (tstar-time1)
    r = r20-r10

    A = norm(v)^2 - a0^2
    B = sum(r .* v) - a0 * a0*tstar
    C = norm(r)^2 - (a*tstar)^2

    if (B > 0 & A >=0) or ((B^2-A*C) < 0)
        return Inf
    else
        return (-B-sqrt(B^2-A*C))/A
    end
end

function jump1(state, k)
  #jump across edge. Given state and index of edge, returns new state
  if k == 1 or k == 3
      ind = 1 #x coordinates
  else
      ind = 2 #y coordinates
  end
  if k == 1 or k ==2
      newval = 0 #right or top edge; jump to left or bottom
  else
      newval = 400 #jump to right or top
  end
  newstate = copy(state)
  newstate[1][ind] = newval
  return newstate
end

function jump2(state1, state2)
  #process collision between sphere 1 and 2

  return state1, state2
end

function growing_disks(n; a0 = 3.2, xlo=0, xhi=400, ylo=0, yhi=400)
  #a(t) = a0*t #DIAMETER GROWTH RATE
  #h = a'(t) = a0
  #we store time, state and partner separately rather than making an actual 'event array'
  current_time = 0
  New = ones(N); Old = twos(N) #initialize 'current index' for each disk
  partner = zeros(N,2); time = zeros(N, 2) #row for disk, column for old vs new
  for i = 1:n
    time[i, 1] = 0 #initiate all times at zero
    partner[i, 1] = NaN #no new partners
    r = rand(2,1).*[xhi-xlo, yhi-ylo]+[xlo, ylo] #random positioning
    v = (rand(2,1) .- 0.5)./sqrt(0.5) #initial velocities between -1 and 1
    state[i, 1] = [ r', v' ]

    time[i,2]=copy(time[i,1]); partner[i,2]=copy(partner[i,1]); state[i,2]=copy(state[i,1])
  end

  while current_time <= end_time
    #current_time, indi = findmin(time[:, New[i]) #nearest new event; ind is index
    #current_time = minimum( vcat(time[:,1][ New .== 1 ], time[:,2][ New .== 2 ]))
    newtimes = zeros(N)
    newtimes[New .== 1] = time[:,1][ New .== 1 ]
    newtimes[New .== 2] = time[:,2][ New .== 2 ]
    current_time, indi = findmin(newtimes)
    New[indi] = Old[indi]
    Old[indi] = 3 - New[indi] #swap new and old
    indi_bar = vcat(1:indi-1, indi+1:N) #indices not equal to i
    A_indi = indi_bar[time[indi_bar, New[indi_bar]] .>= P[indi,:]]
    #P is next interaction time of two spheres
    Pmin = minimum(P[indi, A_indi]) #j = A(ind), A(ind) = {j | 1 <= j <= N, j !=ind, time[j, new[j]] >= P[]ind, j]}
    indj = A_indi[Pmin[1]] #get j value corresponding to this minimum
    Qmin = minimum(Q[indi, :])
    indk = Qmin[1]
    R = min(Pmin, Qmin)

    if R < Inf
      state1 = advance(state[indi, Old[indi]], time[indi, Old[indi]], R, a0)
      if Q < P
        state[indi, New[indi]] = jump1(state1, indk)
        partner[indi, New[indi]] = NaN
      else
        time[indj, New[indj]] = R
        state2 = advance(state[indj, Old[indj]], time[indj, Old[indj]], R, a0)
        state[indi, New[indi]], state[indj, New[indj]] = jump2(state1, state2)
        indm = partner[indj, New[indj]]
        partner[indi, New[indi]] = indj
        partner[indj, New[indj]] = indi
        if !isnan(m) & m != indi #j has collided with i, is no longer ms next partner
          state[indm, New[indm]] = advance(state[indm, New[indm]],
                                          time[indm, Old[indm]],
                                          time[indm, New[indm]])
          partner[indm, New[indm]] = NaN
        end
      end
    end
  end
end

function max_growing_disks(;thresh=10)
  #thresh is min disk ratio; see how many points we can add
  #while still having radius > thresh
end
