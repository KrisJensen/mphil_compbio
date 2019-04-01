function update_2layerold(state, ws, vs, fs, D, cs, A, ptm, ktm; d=dnorm,
                        alpha=alpha, lambda=lambda, gamma=gamma, delta=delta, beta=beta,
                        l=l, mc=mc, mp=mp, muc=muc, mup=mup)

    n = 162
    #ASE
    k = get_state(state)
    zs = 1 ./ (1 .+ exp.(-D[:,k])) #vector
    P = ws[k] + fs' * zs
    P = 1 / (1 + exp(-P)) #probability of right
    #println(P)
    if rand() < P
        yt = 1
    else
        yt = -1
    end
    #ACE
    qs = 1 ./ (1 .+ exp.(-A[:,k])) #vector
    #if ktm == 0 qst1t1 = zeros(n) else qst1t1 = 1 ./ (1 .+ exp.(-A[:,ktm])) end
    pt = vs[k] + cs' * qs
    #if ktm == 0 pt1t1 = 0 else pt1t1 = vs[ktm] + cs' * qst1t1 end

    if state == [0;0;0;0]
        rhatt = 0
    else
        rhatt = gamma*pt - ptm #calculate expected reward at time t
        #rhatt = gamma*pt - pt1t1
    end

    ws[k] += alpha*rhatt*(max(yt,0) - P)
    fs += alpha * rhatt * (max(yt,0) - P)*zs
    D[:,k] += ((0.2*alpha) * rhatt * zs) .* (1 .- zs) .* sign.(fs) * (max(yt,0) - P)

    vs[k] += beta*rhatt
    A[:,k] += ((0.25*beta)*rhatt*qs) .* (1 .- qs) .* sign.(cs)
    #A[:,k] += ((0.25*beta)*rhatt*qst1t1) .* (1 .- qst1t1) .* sign.(cs)
    cs += beta * rhatt * qs
    #cs += beta * rhatt * qst1t1


    state = update_cart(state, yt, l=l, mc=mc, mp=mp, muc=muc, mup=mup) #update states to t+1
    return state, ws, vs, fs, D, cs, A, pt, k, yt, P, zs, qs
end
