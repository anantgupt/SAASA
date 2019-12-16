# Sparse 2D Antenna array design
This repo provides the code for desigining Large Effective Aperture arrays using sparse array of subarrays architecture. Our approach uses multi-objective optimization of key beam parameters to design the array. We provide codes for design and analysis of our algorithm including
1. Synthesizing the array.
2. Performing estimation theoretic bound on the Direction of Arrival Estimation (DoA) performance.
3. Evaluating performance of super-resolution and compressive estimation algorithms.

The steps to use the code are as follows,
## Array Synthesis
1. Run combs9final.m to generate the prefix tree dictionary which is stored in combs9A.mat.
  a. We generated 4 types with different grid sizes, discretization levels.

2. Run case_combs3.m to traverse the tree and evaluate beam parameters for full configurations (leaf node contains config with all Nsub subarrays) which are stored in data_combs9.mat.

3. Analyze the beamparameters and do placement optimization using Lagrangian_optimization.m which returns configuration (subarray centers and poses(up/down)) in <date>.mat.
  a. This script also does local refinement in position and phase to improve cost.
  
![Array Synthesis Demo](demo/Refinements.gif)

## Array analysis
4. Run script_ZZB2.m to analyze the estimation theoretic bounds for array configurations in <date>.mat.
  a. Figures saved in Results2 folder.

## Array performance
5. Run Compressive_Estimation3.m to evaluate performance of algorithms with both compressive and full measurements.
  a. Also evaluates the min DoA separation required for detecting 2 close targets.

# Related Publications
1. Anant Gupta, Upamanyu Madhow, Amin Arbabian, Ali Sadri, ["Design of Large Effective Apertures for Millimeter Wave Systems using a Sparse Array of Subarrays"](https://wcsl.ece.ucsb.edu/sites/default/files/publications/gupta2019design_0.pdf), IEEE Transactions on Signal Processing, vol. 67, no. 24, pp. 6483-6497.
2. A. Gupta, U. Madhow, A. Arbabian and A. Sadri, ["On beam design for sparse arrays of subarrays using multi-objective optimization and estimation-theoretic criteria."](https://wcsl.ece.ucsb.edu/sites/default/files/publications/asilomar17_final.pdf), 51st Asilomar Conference on Signals, Systems and Computers, Nov. 2017, Pacifc Grove, USA.
