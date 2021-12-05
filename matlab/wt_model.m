%% initialize cobra toolbox

updateCobraToolbox
initCobraToolbox

%% import model

% model = readCbModel('../models/mat/iML1515.mat');
model = readCbModel('../models/mat/test.mat');

% outmodel = writeCbModel(model, 'xls', '../results/core_model.xls');

%% genes

genes = cellstr(model.genes);
printLabeledData(model.genes, geneIDs)

% findGeneIDs(model, 'b0116')

% [results ListResults]=findRxnsFromGenes(model, 'b0116', 0, 1)

%% reactions

reactions = cellstr(model.rxns);
rxnIDs = findRxnIDs(model, model.rxns);
printLabeledData(model.rxns, rxnIDs)
r_names = horzcat(reactions, model.rxnNames(rxnIDs));
r_names = horzcat(r_names, model.subSystems);
% r_names = horzcat(r_names, cellstr(num2str(rxnIDs)));

writecell(r_names, '../results/w_reactions.csv')

rxnIDs = findRxnIDs(model, 'glucose')



printRxnFormula(model, 'EX_glyc_e')

[geneList]=findGenesFromRxns(model, 'ENO');
geneList{1:1}

surfNet(model, 'EX_glyc_e')

%% metabolites

metabolites = cellstr(model.mets);
metIDs = findMetIDs(model, model.mets);
printLabeledData(model.mets,metIDs)

% metIDs = findMetIDs(model, 'akg[c]')

% model.metNames(metIDs)

% model.metFormulas(metIDs)

% surfNet(model, 'atp_c')

% rxnIDs = findRxnIDs(model, 'PYK');
% model.grRules(rxnIDs)

% rxnIDs = findRxnIDs(model, 'ATPS4r');
% model.grRules(rxnIDs)

%% constraints

printConstraints(model, -100, +100)

% d-glucose uptake
rxnIDs = findRxnIDs(model, 'EX_glc__D_e');
model.lb(rxnIDs) %-10
model.ub(rxnIDs) %1000

model = changeRxnBounds(model, 'EX_glc__D_e', -10, 'l');
printConstraints(model);

% oxygen uptake
rxnIDs = findRxnIDs(model, 'EX_o2_e');
model.lb(rxnIDs) %-1000
model.ub(rxnIDs) %1000

model = changeRxnBounds(model, 'EX_o2_e', -1000, 'l');
printConstraints(model);



%% objective function

checkObjective(model)

% maximize biomass
printRxnFormula(model, 'BIOMASS_Ec_iJO1366_core_53p95M')


%% fba

model = changeObjective(model,'DM_accoa_c'); % Set the objective function

solverOK = changeCobraSolver('glpk');

FBAsolution = optimizeCbModel(model, 'max') % FBA analysis

