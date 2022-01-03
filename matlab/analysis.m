%% initialize cobra toolbox

initCobraToolbox
% solverOK = changeCobraSolver('glpk');

%% import model

% model = readCbModel('../models/mat/iML1515.mat');
model = readCbModel('../models/mat/GM_iML1515.mat');
% model2 = readCbModel('../models/mat/test2.mat');

verifyModel(model, 'massBalance', true)

%% check model

checkObjective(model)

model = changeRxnBounds(model, 'EX_phleth_e', -14, 'l');
model = changeRxnBounds(model, 'EX_o2_e', -20, 'l');
model = changeRxnBounds(model, 'EX_glc__D_e', -10, 'l');
% model = changeRxnBounds(model, 'EX_so4_e', -20, 'l');


model = changeObjective(model,'DM_phb_c');

%% robustness analysis

figure(1)
% controlRxn = char('DM_accoa_c');
% controlRxn = char('EX_glc__D_e');
controlRxn = char('EX_phleth_e');
% controlRxn = char('EX_o2_e');
[controlFlux, objFlux] = robustnessAnalysis(model, controlRxn)
% xlim([-14 0])
% ylim([0 10])

csvwrite('controlFlux.csv', controlFlux)
csvwrite('objFlux.csv', objFlux)

%% double robustness analysis

figure(1)
controlRxn1 = char('EX_glc__D_e');
controlRxn2 = char('EX_phleth_e');
doubleRobustnessAnalysis(model, controlRxn1, controlRxn2)

%%
figure(2)
controlRxn1 = char('EX_phleth_e');
controlRxn2 = char('EX_o2_e');
[controlFlux1, controlFlux2, objFlux] = doubleRobustnessAnalysis(model, controlRxn1, controlRxn2);

%% optknock for acetyl coa maximization (currently broken)

biomass = 'BIOMASS_Ec_iML1515_core_75p37M';

threshold = 5;

selectedRxnList = {'GLCabcpp'; 'GLCptspp'; 'HEX1'; 'PGI'; 'PFK'; 'FBA'; 'TPI'; 'GAPD';
'PGK'; 'PGM'; 'ENO'; 'PYK'; 'LDH_D'; 'PFL'; 'ALCD2x'; 'PTAr'; 'ACKr';
'G6PDH2r'; 'PGL'; 'GND'; 'RPI'; 'RPE'; 'TKT1'; 'TALA'; 'TKT2'; 'FUM';
'FRD2'; 'SUCOAS'; 'AKGDH'; 'ACONTa'; 'ACONTb'; 'ICDHyr'; 'CS'; 'MDH';
'MDH2'; 'MDH3'; 'ACALD'};

solution = optimizeCbModel(model, 'max');

growth_rate = solution.f;
accoa_yield = solution.x(strcmp(model.rxns, 'DM_accoa_c'));
phb_yield = solution.x(strcmp(model.rxns, 'DM_phb_c'));

fprintf('growth rate prior to optknock: %.1f \n', growth_rate);
fprintf('accoa yield prior to optknock: %.1f \n', accoa_yield);
fprintf('phb yield prior to optknock: %.1f \n', phb_yield);

fprintf('\n...finding optKnock sets of size 2 or less...\n\n')

% Set optKnock options
% The exchange of succinate will be the objective of the outer problem
options = struct('targetRxn', 'EX_succ_e', 'numDel', 2);
% We will impose that biomass be at least 50% of the biomass of wild-type
constrOpt = struct('rxnList', {{biomass}},'values', 0.5*solution.f, 'sense', 'G');
% We will try to find 10 optKnock sets of a maximun length of 2
previousSolutions = cell(10, 1);
contPreviousSolutions = 1;
nIter = 1;
while nIter < threshold
fprintf('...performing optKnock analysis...\n')
if isempty(previousSolutions{1})
optKnockSol = OptKnock(model, selectedRxnList, options, constrOpt);
else
optKnockSol = OptKnock(model, selectedRxnList, options, constrOpt, previousSolutions, 1);
end
% determine succinate production and growth rate after optimization
growth_rateM1 = optKnockSol.fluxes(strcmp(model.rxns, biomass));
accoa_yieldM1 = optKnockSol.fluxes(strcmp(model.rxns, 'DM_accoa_c'));
phb_yieldM1 = optKnockSol.fluxes(strcmp(model.rxns, 'DM_phb_c'));
setM1 = optKnockSol.rxnList;
if ~isempty(setM1)
previousSolutions{contPreviousSolutions} = setM1;
contPreviousSolutions = contPreviousSolutions + 1;
%printing results
fprintf('optKnock found a optKnock set of large %d composed by ', length(setM1));
for j = 1:length(setM1)
if j == 1
fprintf('%s', setM1{j});
elseif j == length(setM1)
fprintf(' and %s', setM1{j});
else
fprintf(', %s', setM1{j});
end
end
fprintf('\n');
fprintf('growth rate after optimization: %.2f \n', growth_rateM1);
fprintf('accoa yield after optimization: %.2f \n', accoa_yieldM1);
fprintf('phb yield after optimization: %.2f \n', phb_yieldM1);
fprintf('...Performing coupling analysis...\n');
[type, maxGrowth, maxProd, minProd] = analyzeOptKnock(model, setM1, 'DM_accoa_c');
fprintf('solution is of type: %s\n', type);
fprintf('maximun growth rate given the optKnock set: %.2f\n', maxGrowth);
fprintf(['maximun and minimun production of acetyl-coa given the optKnock set: '
'%.2f and %.2f, respectively \n\n'], minProd, maxProd);
if strcmp(type, 'growth coupled')
singleProductionEnvelope(model, setM1, 'DM_accoa_c', biomass, 'savePlot', 1, 'showPlot', 1, ...
'fileName', ['accoa_dm1_' num2str(nIter)], 'outputFolder', 'OptKnockResults');
end,
else
if nIter == 1
fprintf('optKnock was not able to found an optKnock set\n');
else
fprintf('optKnock was not able to found additional optKnock sets\n');
end
break;
end
nIter = nIter + 1;
end


