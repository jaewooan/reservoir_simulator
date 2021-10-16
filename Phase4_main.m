close all;
clear;
clc;
runTime=0;
%% Phase 1 DFM
% filename = 'ref_dpdk3.yaml';
% [fluid_, grid_, connection_, equation_, well_, output_, frac_, runTime] = ...
%     simulator(filename);
filename = 'ref_edfm.yaml';
[fluid_, grid_, connection_, equation_, well_, output_, frac_, runTime] = ...
    simulator(filename);
% make new input file: specifically dfeine fracture information
% modify initialize
% modify connection list

% 
%% Phase 2 EDFM
filename2 = 'ref_dfm_edfm.yaml';
[fluid2_, grid2_, connection2_, equation2_, well2_, output2_, frac2_, runTime2] = ...
    simulator(filename2);
% % 
%% Check it is fit to the Eclipse
filename = 'A/PHASE4A';
iFigure = plotResult(2, output_, filename, fluid_, grid_, well_, frac_,...
                      output2_, filename2, fluid2_, grid2_, well2_, frac2_);
% % %                     
% %% Phase 1 DFM
% filename = 'ref_dfm_simple_c.yaml';
% [fluid_, grid_, connection_, equation_, well_, output_, frac_, runTime] = ...
%     simulator(filename);
% % make new input file: specifically dfeine fracture information
% % modify initialize
% % modify connection list
% 
% 
% %% Phase 2 EDFM
% filename2 = 'ref_edfm_simple_c.yaml';
% [fluid2_, grid2_, connection2_, equation2_, well2_, output2_, frac2_, runTime2] = ...
%     simulator(filename2);
% 
% %% Check it is fit to the Eclipse
% filename = 'A/PHASE4A';
% iFigure = plotResult(iFigure, output_, filename, fluid_, grid_, well_, frac_,...
%                         output2_, filename2, fluid2_, grid2_, well2_, frac2_);