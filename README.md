# 714-Final-Project

# Spectral Splitting Methods
spectral_Lie_KdV2.m 

spectral_Strang_KdV3.m

spectral_SWSS_KdV.m 

All of these take in the arguments a,b,N,T where [a,b] is the interval, N is the number of grid points, and T is the time up to which we solve. They use the different splitting methods all in Fourier space. Initial data must be changed within each file. 

# Other Splitting Methods
KdV_Lie_LF.m

KdV_Lie_Godunov.m

Both of these take in the arguments a,b,N,T where [a,b] is the interval, N is the number of grid points, and T is the time up to which we solve. They use the Lie splitting but with Burgers' equation solved in real space. Initial data must be changed within each file. 

F.m - helper function which calculated Godunov flux

# Solution Quantities
mass.m

momentum.m

Both of these take in a solution u and output an array which length is the number of time steps and calculates the mass or momentum of the solution at each time step. 

# Post Processing
convergence.m - Calculates errors of solution to an exact solution. Currently set for a single soliton solution. Also calculates M_diff and p_diff values.

multiple_plotter.m - Plots solution at various times. 
