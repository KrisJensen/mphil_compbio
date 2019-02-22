
using PyCall, PyPlot, LinearAlgebra, Random, DelimitedFiles, LaTeXStrings
#use the matplotlib 'quiver' function for gradient fields



MEE = 1.25
MIE = 1
MII = -1
MEI = -1
gE = -10
gI = 10
tE = 10
dt = 0.01


function simulate(;tI=70, dt=0.1, T=10000, v0s = [15,10])
    N = Int(T/dt)+1
    vs = zeros(N, 2)
    vs[1,:] = v0s

    for i in 2:N
        vE, vI = vs[i-1,:]
        vs[i,1] = vE + dt/tE*( max((MEE*vE + MEI * vI - gE),0) - vE )
        vs[i,2] = vI + dt/tI*( max((MII*vI + MIE * vE - gI), 0) - vI )
    end
    return vs
end

function plotsim(vs; fname="trajectcon.png", tI=70)
    figure(figsize=(5,3.5))
    plot(vs[:,1], vs[:,2], "k-")
    vEmax = 134#maximum(vs[:,1])*1.05
    vImax = 52#maximum(vs[:,2])*1.05

    xlim([0,vEmax])
    ylim([0,vImax])
    vEs = 1:7:vEmax
    vIs = 1:3:vImax
    nvI = [get_vI_null(vE) for vE in vEs]
    nvE = [get_vE_null(vE) for vE in vEs]
    plot(vEs, nvI, "r--")
    plot(vEs, nvE, "b--")

    U = zeros(length(vIs), length(vEs))
    V = zeros(length(vIs), length(vEs))
    X = zeros(length(vIs), length(vEs))
    Y = zeros(length(vIs), length(vEs))

    for i in 1:length(vIs)
        for j in 1:length(vEs)
            vI=vIs[i]
            vE = vEs[j]
            U[i,j] = 1/tE*( max((MEE*vE + MEI * vI - gE),0) - vE )
            V[i,j] = 1/tI*( max((MII*vI + MIE * vE - gI), 0) - vI )
            X[i,j] = vE
            Y[i,j] = vI
        end
    end
    quiver(X, Y, U, V, pivot="tail")

    xlabel(L"v_E (Hz)")
    ylabel(L"v_I (Hz)")
    savefig(fname, bbox_inches="tight", dpi=300)
    show()
    close()
end

vs = simulate(tI=75, v0s=[20,10])
plotsim(vs, tI=75)
vsdiv = simulate(tI=85, v0s=[60,21])
plotsim(vsdiv, fname="trajectdiv.png", tI=85)
vsn = simulate(tI=80)
plotsim(vsn, fname="trajectn.png", tI=80)

#dvE = dt/tE*( (MEE*vE + MEI * vI - gE) - vE )
#dvI = dt/tI*( (MII*vI + MIE * vE - gI) - vI )




function get_vE_null(vE)
    return (gE+vE*(1-MEE))/MEI
end
function get_vI_null(vE)
    return (MIE*vE-gI)/(1-MII)
end
function plot_nulls()
    vEs = 1:0.1:100
    nvI = [get_vI_null(vE) for vE in vEs]
    nvE = [get_vE_null(vE) for vE in vEs]
    figure(figsize=(5,3.5))
    plot(vEs, nvI, "r--")
    plot(vEs, nvE, "b--")
    xlim([0,maximum(vEs)])
    ylim([0,maximum(nvI)])
    xlabel(L"\tau_E (Hz)")
    ylabel(L"\tau_I (Hz)")
    show()
    savefig("nullclines.png", bbox_inches="tight", dpi=300)
    close()
end
#plot_nulls()

function get_lambds(tI)

    D = ( (MEE-1)/tE + (MII-1)/tI )^2 + 4*MEI*MIE/(tE*tI)

    if D >= 0
        lambd1 = 0.5*( (MEE-1)/tE + (MII-1)/tI + sqrt(D) )
        lambd2 = 0.5*( (MEE-1)/tE + (MII-1)/tI - sqrt(D) )
        return([lambd1, lambd2], [0,0])
    else
        lambd1 = 0.5*( (MEE-1)/tE + (MII-1)/tI )
        lambd2 = 0.5*( (MEE-1)/tE + (MII-1)/tI )
        return([lambd1, lambd2], [0.5*sqrt(-D), -0.5*sqrt(-D)])
    end

end



function test_tIs()
    tIs = 0.1:0.1:200
    lambd1s = [get_lambds(tI)[1][1] for tI in tIs]
    lambd2s = [get_lambds(tI)[1][2] for tI in tIs]
    figure(figsize=(5,3.5))
    plot(tIs, lambd1s, "r--")
    plot(tIs, lambd2s, "b--")
    plot([0, tIs[end]], [0, 0], "k:")
    plot([80, 80], [-2, 2], "k:")
    ylim([-0.12, 0.025])
    xlabel(L"\tau_I (s)")
    ylabel(L"\lambda (Hz)")
    legend([L"\lambda_1", L"\lambda_2"])
    savefig("ti_scan.png", bbox_inches="tight", dpi=300)
    close()
    print(maximum(lambd1s), maximum(lambd2s))
end
#test_tIs()

vIfp = (gE/(MEE-1)-gI/MIE) / (MEI/(MEE-1)-(MII-1)/MIE)
vEfp = gE/(MEE-1) - MEI/(MEE-1)*vIfp
