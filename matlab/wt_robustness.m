%% initialize cobra toolbox

initCobraToolbox
solverOK = changeCobraSolver('glpk');

%% import model

model = readCbModel('../models/mat/iML1515.mat');
% model = readCbModel('../models/mat/test.mat');

%% check model

checkObjective(model)

%% robustness

% controlRxn = char('DM_accoa_c');
controlRxn = char('EX_ac_e');
robustnessAnalysis(model, 'EX_ac_e')
