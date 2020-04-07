%This script creates the variables through which the required parameters 
%and files are inputted to the metagenomic pipeline (MgPipe). Input 
%variables should be changed by the user according to what specified in the 
%documentation. Running this script will automatically launch the pipeline. 

% Federico Baldini, 2017-2018

initCobraToolbox()
%REQUIRED INPUT VARIABLES

% path to microbe models
modPath='YOUR_PATH_TO_AGORA\';
% path where to save results
resPath='YOUR PATH TO RESULT FOLDER\';
% path to where the COBRA Toolbox is located
global CBTDIR
toolboxPath=CBTDIR;
% path to and name of the file with dietary information.
dietFilePath=[CBTDIR filesep 'papers' filesep '2018_microbiomeModelingToolbox' filesep 'resources' filesep 'AverageEuropeanDiet'];
% path to and name of the file with abundance information.
abunFilePath=[CBTDIR filesep 'papers' filesep '2018_microbiomeModelingToolbox' filesep 'examples' filesep 'normCoverage.csv'];
% path to csv file for stratification criteria (if empty or not existent no criteria is used)
indInfoFilePath='none'
% name of objective function of organisms 
objre={'EX_biomass(e)'};
%the output is vectorized picture, change to '-dpng' for .png
figForm = '-depsc'
% number of cores dedicated for parallelization 
numWorkers = 3;
% autofix for names mismatch
autoFix = true;
% if outputs in open formats should be produced for each section 
compMod = false; 
% to enable also rich diet simulations 
rDiet =false;
% if if to use an external solver and save models with diet
extSolve = false;
% the type of FVA function to use to solve
fvaType=true;
% To tourn off the autorun to be able to manually execute each part of the pipeline.
autorun=false;
%END OF REQUIRED INPUT VARIABLES

%%
%PIPELINE LAUNCHER 
[init,modPath,toolboxPath,resPath,dietFilePath,abunFilePath,indInfoFilePath,objre,figForm,numWorkers,autoFix,compMod,rDiet,extSolve,fvaType,autorun]= initMgPipe(modPath, toolboxPath, resPath, dietFilePath, abunFilePath,indInfoFilePath, objre, figForm, numWorkers, autoFix, compMod, rDiet,extSolve,fvaType,autorun);


