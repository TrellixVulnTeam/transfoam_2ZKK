%% initialize cobra toolbox

initCobraToolbox
% changeCobraSolver('glpk');
changeCobraSolver('gurobi6');

%% import model

cobra_model = readCbModel('../models/mat/GM_iML1515.mat');

cobra_model = changeRxnBounds(cobra_model, 'EX_phleth_e', -10, 'l');
cobra_model = changeRxnBounds(cobra_model, 'EX_o2_e', -20, 'l');
cobra_model = changeRxnBounds(cobra_model, 'EX_glc__D_e', 0, 'l');

% cobra_model = changeObjective(cobra_model,'DM_phb_c');

%% optpipe

%Disable warnings
warning('off', 'all');
%TODO: warnings should be handled

% parameter settings
BM_rxn='BIOMASS_Ec_iML1515_core_75p37M';
target_rxn='DM_phb_c';
max_KOs=3;
branch='';

[optknock_results,robokod_results,optgene_results,screening_results_max,screening_results_min] = ...
optPipe(cobra_model,BM_rxn,target_rxn,max_KOs,branch)

success = true;

%Re-enable warnings
warning('on', 'all');