#Primary file for simulations used in spa3. Calls functions implemented
#In remaining files to run simulations.

using Random, Distributions, PyCall, PyPlot, LinearAlgebra, DelimitedFiles, Plots
@pyimport matplotlib.patches as patch

#We start by loading functions from auxillary files
include("plots.jl") #functions for plotting everything
include("dmin2d.jl") #main dmin2d function
include("RI50.jl") #functions for calculating regularity indices and u values
include("saturation.jl") #functions for quantifying how many points we can add to a grid
include("LS_model.jl") #functions used in Lubachevsky-Stillinger simulations

####Now run actual code

#generate and plot grids for part I
for i in 1:2 points = dmin2d(200, 30, 5, 200, 1000, 100, 900, Print=false, Plot=true) end

RI = calc_RI(points) #test calculation of regularity index
run_sims(;N = 1000, n=200, m=0, s=0, xlo=0, xhi=1000, ylo=0, yhi=1000) #calculate RI50
results = vary_params() #investigate effect of n and shape

ref = readdlm("spa3_real.dat", Float64) #load reference points
calc_sim(10, 5, ref) #test similarity calculation
result_mat = scan_ms(ref=ref) #map u1 as a function of m and s
results = steepest_descent() #run steepest descent for optimum (m, s)
sd = test_sd() #find standard deviation of u1 for 10 simulations
plot_descent() #plot result of steepest descent and project onto u1 landscape
get_opt_means(;sds = 0:0.2:5.4) #find optimum m for for each s to minimize u1

plot_packings() #plot square and hexagonal packings
ns = repeat_max() #perform 10 bisection searches to find maximum dmin2d packing
ns = repeat_max(method="bd") #perform 10 bisection searches to find maximum birth-death packing

LS_simulation(;n=12, a0 = 2.5, xlo=0, xhi=420, ylo=0, yhi=420) #pack 12 spheres
LS_simulation(;n=24, a0 = 2.5, xlo=0, xhi=420, ylo=0, yhi=420) #pack 24 spheres
LS_simulation(;n=24, a0 = 4.5, xlo=0, xhi=420, ylo=0, yhi=420) #pack 24 spheres with high a0
repeat_max_LS() #find maximum packing with Lubachevsky-Stillinger model
