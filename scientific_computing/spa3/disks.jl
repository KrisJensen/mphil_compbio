function find_next_event(state, xlo, xhi, ylo, yhi)

end

function advance(state0, t0, t1, a0)
  #compute the state after time t ignoring collisions
  #update velocity and position
  r = state[1], v=state[2]

  return state1

end

function jump1(state, indk)
  #jump across edge

  return state
end

function jump2(state1, state2)
  #process collision between sphere 1 and 2

  return state1, state2
end

function growing_disks(n; a0 = 3.2, xlo=0, xhi=400, ylo=0, yhi=400)
  #a(t) = a0*t
  #h = a'(t) = a0
  current_time = 0
  New = ones(N); Old = twos(N) #initialize 'current index' for each disk
  partner = zeros(N,2); time = zeros(N, 2) #row for disk, column for old vs new
  for i = 1:n
    time[i, 1] = 0 #initiate all times at zero
    partner[i, 1] = NaN #no new partners
    r = rand(2,1).*[xhi-xlo, yhi-ylo]+[xlo, ylo] #random positioning
    v = (rand(2,1) .- 0.5)./sqrt(0.5) #initial velocities between -1 and 1
    state[i, 1] = [ r', v' ]

    event[i,2] = event[i,1] #what is this???
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
