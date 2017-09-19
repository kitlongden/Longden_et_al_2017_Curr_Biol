%% To plot local motion receptive fields of VT1, VT2, VT3 and Hu cells, 
%% run <Plot_local_motion_receptive_field.m> in this folder. It will load the data 
%% from the ../Data folder. The function takes two arguments, the recording
%% identification number <rec_id> and the cell type <cell_type>. 
%% Cell type 1 = VT1 wide field maps, 
%%           2 = VT1 ipsilateral ventral receptive field maps
%%           3 = VT2
%%           4 = VT3
%%           5 = Hu.
%% 
%% To calculate the VT1 ipsilateral fits, as used in Figure 1, 
%% run <Calculate_VT1_ipsilateral_optic_flow_fits.m>. 
%% This calls the <Fit_rot_trans_ipsi.m>, function which calculates the fit.
%% In doing so, the functions <intermap.map> is called to interpolate the raw data
%% maps, and <Generate_rotation_flow_map.m> and <Generate_rotation_flow_map.m> are
%% called to generate the rotation and translation optic flow maps.
%% 
%% To calculate the optic flow field fields using both ipsilateral and contralateral
%% visual fields, run <Calculate_VT2_VT3_Hu_contralateral_optic_flow_fits.m>.
%%
%% To calculate the log ISI distribution of VT1, run <Calculate_ISI_for_Figure_1.m>.
%% 
%% Kit Longden
%% 6th September 2017
