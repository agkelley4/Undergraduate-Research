# Flow of Deformable Cells Undergraduate Research

Spherical Harmonic Equation Sheet show the equations in spherical coordinates used to model spherical harmonics. Phi is from 0 to pi and theta is from 0 to 2pi.

Spherical Harmonic Plot table shows the basic spherical harmonics with coefficient 1

Run plot_sph_harm.py with coeffic_to_plot.py to model spherical harmonic linear combinations. Format CSV file to el,m,alpha \n "el value","m value","coefficient value"

Volume Integration Code.nb is Mathematica code used to take the volume integral of spherical harmonics and linear combinations of spherical harmonics. 

Used linear_comb.m to create the figure, constraint_Y00_Y10_with_sph_harm.png. Shows all combinations of alpha and beta coefficients that will hold the following condition, alphaY00 + betaY10 = Y00. In order to calculate the volume for alphaY10 + betaY10 we had to 

Used Monte Carlo integration to verify our calculations for volume of alphaY00 + betaY10 = Y00. Found under montecarloy00_cartesian_jun6th2024.m. Inorder to verify if each randomized point was inside the shape, we first converted the cartesian point to spherical point. Next we calculated the radius value of our shape using the same theta and phi value from our randomized point. This point gives us the shell of the shape. Knowing all spherical harmonics lie inside the origin, we then drew a line from the shell to the origin. Finally we can check if our inital randomized point is inside the shape by checking if its on the line or not.


I used a parameter sweep in COMSOL to allocate data for forces and torque in XYZ directions for every combination of phi and theta. plt_Fp_final.py and plt_Tp_final.py creates a plot showing the forces and torques, respectivly, compared to each combination of phi and theta.
Use force_ellipsoid_theta_psi_fullsweep.txt with plt_Fp_final.py and plt_Tp_final.py.
Output example is found at Example_plt_Fp_final.png

Using another parameter sweep in COMSOL, I found the values of forces and torque with respect to every combination of velocity and angular velocity of the cell. plt_Ompy_final.py and plt_Up.py use this data to create plots showing how forces and torque compare to Ompy and Up
Use force_ellipsoid_Ompy_Up_fullsweep.txt with plt_Ompy_final.py and plt_Up.py.
Output example is found at Example_plt_Up.png

Using information gathered from parameter sweeps, optimize_Fz_T.m uses inputs of phi, theta, Up,
