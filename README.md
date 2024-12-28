# Flow of Deformable Cells Undergraduate Research

/Spherical Harmonic Plotter contains plots of deformable cells and code used to create plots. Contains table of first nine spherical harmonics.

/Constant Volume Linear Combinations contains solved plots of linear combinations of spherical harmonics that constant volume with the first spherical harmonic Y00. Used Monte Carlo integration to verify our calculations for volume of alphaY00 + betaY10 = Y00. Found under montecarloy00_cartesian_jun6th2024.m. In order to verify if each randomized point was inside the shape, we first converted the cartesian point to spherical point. Next we calculated the radius value of our shape using the same theta and phi value from our randomized point. This point gives us the shell of the shape. Knowing all spherical harmonics lie inside the origin, we then drew a line from the shell to the origin. Finally we can check if our initial randomized point is inside the shape by checking if its on the line or not.

/Flow Optimization contains a data analysis on 

I used a parameter sweep in COMSOL to allocate data for forces and torque in XYZ directions for every combination of phi and theta. plt_Fp_final.py and plt_Tp_final.py creates a plot showing the forces and torques, respectivly, compared to each combination of phi and theta.
Use force_ellipsoid_theta_psi_fullsweep.txt with plt_Fp_final.py and plt_Tp_final.py.
Output example is found at Example_plt_Fp_final.png

Using another parameter sweep in COMSOL, I found the values of forces and torque with respect to every combination of velocity and angular velocity of the cell. plt_Ompy_final.py and plt_Up.py use this data to create plots showing how forces and torque compare to Ompy and Up
Use force_ellipsoid_Ompy_Up_fullsweep.txt with plt_Ompy_final.py and plt_Up.py.
Output example is found at Example_plt_Up.png

Using information gathered from parameter sweeps, optimize_Fz_T.m uses inputs of phi, theta, Up,
