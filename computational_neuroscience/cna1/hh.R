library(deSolve)

cm. = 10 #nF/mm2
V. = -65 #mV
m. = 0.0529 #unitless gating variables
h. = 0.5961
n. = 0.3177
dt. = 0.1 #ms  
Ie_A. = 200 #nA/mm2   (nA/mm2 / nF/mm2 = A/F = V/s))

#physical constants from TN
gL = 3 #uS/mm2 (from TN)
gK = 360 #uS/mm2
gNa = 1200 #uS/mm2 gna*ENa: us*mV/mm2 = nA/mm2
EL = -54.387 #mV
EK = -77 #mV
ENa = 50 #mV

#update

calc_params = function(state){
  V = state['V']; n = state['n']; m = state['m']; h = state['h']
  #potassium
  an = 0.01*(V+55)/(1-exp(-0.1*(V+55)))
  bn = 0.125*exp(-0.0125*(V+65))
  tn = 1 / (an+bn)
  ninf = an * tn
  #sodium
  am = 0.1*(V+40)/(1-exp(-0.1*(V+40)))
  bm = 4*exp(-0.0556*(V+65))
  tm = 1 / (am+bm)
  minf = am * tm
  ah = 0.07*exp(-0.05*(V+65))
  bh = 1/(1+exp(-0.1*(V+35)))
  th = 1 / (ah+bh)
  hinf = ah * th
  #print(V)
  im = gL*(V-EL) + gK*(n^4)*(V-EK) + gNa*(m^3)*h*(V-ENa)
  return(list(tn=tn, tm=tm, th=th, im=im))
}

hh = function(t, state, parameters) {
  V = state['V']; n = state['n']; m = state['m']; h = state['h']
  #potassium
  an = 0.01*(V+55)/(1-exp(-0.1*(V+55)))
  bn = 0.125*exp(-0.0125*(V+65))
  tn = 1 / (an+bn)
  ninf = an * tn
  #sodium
  am = 0.1*(V+40)/(1-exp(-0.1*(V+40)))
  bm = 4*exp(-0.0556*(V+65))
  tm = 1 / (am+bm)
  minf = am * tm
  ah = 0.07*exp(-0.05*(V+65))
  bh = 1/(1+exp(-0.1*(V+35)))
  th = 1 / (ah+bh)
  hinf = ah * th
  #print(V)
  im = gL*(V-EL) + gK*(n^4)*(V-EK) + gNa*(m^3)*h*(V-ENa)
  
  with(as.list(c(state, parameters)),{
    
    # rate of change
    dV = (Ie_A-im)/cm
    dn = (ninf-n)/tn
    dm = (minf-m)/tm
    dh = (hinf-h)/th
    
    # return the rate of change
    list(c(dV, dn, dm, dh))
  })# end with(as.list ...
}

run_hh = function(cm=cm., Ie_A=Ie_A., V=V., n=n., m=m., h=h.,
                  dt=dt., tstop=30, Plot=TRUE, jump='none'){
  
  parameters = c(cm = cm,
              Ie_A = Ie_A)
  state = c(V = V,
          n = n,
          m = m,
          h = h)

  if (jump == 'none'){
    times <- seq(0, tstop, by = dt)
    out <- ode(y = state, times = times, func = hh, parms = parameters)
  }else{
    times = seq(0, jump[1], by=dt)
    out1 <- ode(y = state, times = times, func = hh, parms = parameters)
    times = seq(jump[1], tstop, by=dt)
    parameters['Ie_A'] = jump[2]
    state[1:4] = out1[dim(out1)[1], 2:5]
    out2 <- ode(y = state, times = times, func = hh, parms = parameters)
    out = rbind(out1,out2)
  }
  if(Plot){
    par(oma = c(0, 0, 3, 0))
    plot(out, xlab = "time", ylab = "-")
    mtext(outer = TRUE, side = 3, "Hodgkin Huxley Model", cex = 1.5)
    par(mfrow=c(1,1), oma = c(1,1,1,1))
  }
  return(out)
}


get_rate = function(times, Vs){
  spikes = c()
  for (i in 2:(length(Vs)-1)){
    if ( Vs[i] > Vs[i-1] & Vs[i]>Vs[i+1] & Vs[i] > 0){spikes = c(spikes, times[i])}
  }
  rate = 0
  if(length(spikes) >= 2){
    int = spikes[length(spikes)]-spikes[length(spikes)-1]
    ratio = length(spikes)*int / times[length(times)]
    if (ratio > 0.5 & ratio < 2){ #need sustained firing
      rate = 1000/(int)
    }
  }
  return( list(rate, spikes) )
}

################investigate rates
scanrates = function(Ie_As = seq(0,500,5), fname="firing_rates.pdf"){
  rates = c()
  for (IeA in Ie_As){
    tstop = 400
    N = tstop/dt
    out = run_hh(Ie_A=IeA, tstop=tstop, Plot=FALSE)
    a = get_rate(out[,'time'], out[,'V'])
    rates = c(rates, a[[1]])
  }
  pdf(fname, height=6, width=10)
  par(mfrow=c(1,1), oma = c(1,1,1,1))
  plot(Ie_As, rates, xlab = "Ie/A (nA/mm2)", ylab="firing rate (hz)", pch=4, cex=0.5)
  dev.off()
}


#################plot parameters
plotstates = function(jump = 'none', fname='hh_', Ie_A=200, tstop=30, lineval = 5.0){
  out = run_hh(Plot=TRUE, Ie_A=Ie_A, jump=jump, tstop=tstop)
  ylabs = c("V (mV)", "n", "m", "h")
  state = out[dim(out)[1], 2:5]
  print(names(state))
  for (i in 1:4){
    name = names(state)[i]
    pdf(paste0(fname, name, '.pdf'), height=6, width=10)
    par(mfrow=c(1,1), oma = c(1,1,1,1))
    plot(out[,'time'], out[, name], type="l", xlab="time (ms)", ylab =  ylabs[i])
    abline(v=lineval, lty=1, lwd=0.4)
    dev.off()
  }
  
  pdf(paste0(fname, 'gating.pdf'), height=6, width=10)
  par(mfrow=c(1,1), oma = c(1,1,1,1))
  plot(out[,'time'], out[, 'n'], type='l', col="red", ylim=c(0,1), xlab="time (ms)", ylab="")
  lines(out[,'time'], out[, 'm'], lty=2, xlab="time (ms)", col="blue")
  lines(out[,'time'], out[, 'h'], lty=3, xlab="time (ms)")
  legend(27.5, 0.95, c('n', 'm', 'h'), lty=c(1,2,3), col=c("red", "blue", "black"))
  abline(v=lineval, lty=1, lwd=0.4) 
  dev.off()
  
  return(out)
}
out = plotstates(lineval = 13.6)


#scanrates()
#scanrates(Ie_As = seq(60, 65, 0.1), fname="firing_rates_finegrained.pdf")

#out1 = plotstates(Ie_A=-50, jump=c(5, 0), fname='hh_jump_', tstop=40, lineval=5.0)
