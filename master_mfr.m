% Runs available functions for plotting of reach and grasp population analyses. 
% Make sure the following are on the path: 
%
% marmo_grasp_model.mat
% marmo_reach_model.mat
% reachsim_preddelay.mat
% marm_regress_1.mat
%
% Shaw,L, Wang KH, Mitchell, J (2023) Fast Prediction in Marmoset Reach-to-Grasp Movements for Dynamic Prey.
%
% Jude Mitchell, Kuan Hong Wang, and Luke Shaw 4/2023
% MATLAB R2022b
%
% Reaching data structure marmo_reach_model.mat available at
% https://doi.org/10.5281/zenodo.7869286

%Figure 2D, 2E, 2F
grasp_cluster_mfr();

%Figure 3C
plot_normalized_velocity_mfr();

%Figure 3G, 3H
pursuit_model_compare_mfr();

%Figure 4A, 4B
estimate_delay_mfr();

%Figure 4E
regression_model_plot_mfr(1);
%regression_model_mfr(n) to run analysis

%Figure 3I, 4C
pursuit_model_pred_delay_mfr();
