function cbmodel= CNAcna2cobra(cnap)
% CellNetAnalyzer API function 'CNAcna2cobra'
%
% Usage:  cbmodel = CNAcna2cobra(cnap)
%
% Input: cnap is a CNA mass-flow project structure
%
% Output: cbmodel is a COBRA Toolbox model.
%
% Creates a COBRA Model from a CNA project.


% create new CNA Project with default marcomolecule proportions
cnap2= CNAgetMFNetwork(cnap);

cbmodel = struct();
cbmodel.description = 'COBRA model converted from a CNA project';
cbmodel.rxns = cellstr(cnap2.reacID);
cbmodel.mets = cellstr(cnap2.specID);
cbmodel.S =  cnap2.stoichMat;
cbmodel.lb =  cnap2.reacMin;
cbmodel.ub =  cnap2.reacMax;
cbmodel.c =   cnap2.objFunc; 
cbmodel.osense = 1;% CNA minimizes!
cbmodel.osenseStr = 'min'; %Set both fields
cbmodel.rxnNames = cellstr(cnap2.reacID);
cbmodel.metNames = cellstr(cnap2.specLongName);

% optional fields of COBRA model
%  cbmodel.subSystems =      'Subsystem name for each reaction (opt)  ';
%  cbmodel.grRules =         'Gene-reaction association rule for each reaction (opt) ';
%  cbmodel.rules =           'Gene-reaction association rule in computable form (opt)';
%  cbmodel.rxnGeneMat =      'Reaction-to-gene mapping in sparse matrix form (opt)   ';
%  cbmodel.genes =           'List of all genes (opt)      ';
%  cbmodel.metFormulas =     'Metabolite chemical formula (opt)';

% find indices of external metabolites
indices = find(cnap2.specExternal);
disp(['Removing ',num2str(numel(indices)),' external metabolites']);

% remove the external metabolite rows from stoichiometry matrix
cbmodel.S(indices,:)=[];

% create sparse matrix
cbmodel.S =  sparse(cbmodel.S);

% remove the external metabolite ids/Names
cbmodel.mets(indices,:)=[];
cbmodel.metNames(indices,:)=[];
cbmodel.b=zeros(size(cbmodel.S,1),1);



