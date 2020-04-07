function [baseMetNames, compSymbols, uniqueMetNames, uniqueCompSymbols] = parseMetNames(metNames)
% Figures out the base metabolite names and compartments for each metabolite
%
% USAGE:
%
%    [baseMetNames, compSymbols, uniqueMetNames, uniqueCompSymbols] = parseMetNames(metNames)
%
% INPUT:
%    metNames:             List of metabolite names
%
% OUTPUTS:
%    baseMetNames:         List of met names without compartment symbol
%    compSymbols:          Compartment symbols for each metabolite
%    uniqueMetNames:       Unique metabolite names (w/o comp symbol)
%    uniqueCompSymbols:    Unique compartment symbols
%
% Metabolite names should describe the compartment assignment in the
% form "metName[compName]" 
%
% .. Author: - Markus Herrgard 10/4/06
%            - Thomas Pfau Speedup and cleanup Oct 2017
if ~iscell(metNames)
    metNames = {metNames};
end

try
    data = cellfun(@(x) regexp(x,'^(?<metNames>.*)\[(?<compSymbols>[^\[*])\]$','names'),metNames);
catch
    %If the above doesn't work, its likely, that we either have an odd
    %compartment symbol (e.g. (), or that we have a non compartmented
    %entry.
    %Lets see if we have a ()compartment symbol
    try
        data = cellfun(@(x) regexp(x,'^(?<metNames>.*)\((?<compSymbols>[^\[*])\)$','names'),metNames);
    catch
        %No, we don't lets assume, we only have a metabolite name
        data = cellfun(@(x) regexp(x,'^(?<metNames>.*)[\(\[]*(?<compSymbols>.*)[^\[\(]*[\]\)]*$','names'),metNames);
    end
end
baseMetNames = columnVector({data.metNames});
compSymbols = columnVector({data.compSymbols});
uniqueCompSymbols = unique(compSymbols);
uniqueMetNames = unique(baseMetNames);
