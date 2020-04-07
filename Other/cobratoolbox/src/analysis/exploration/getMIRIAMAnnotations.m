function annotations = getMIRIAMAnnotations(model, type, varargin)
% Get the MIRIAM Annotations a struct array of annotation structs for the given type
%
% USAGE:
%    annotations = getMIRIAMAnnotations(model,field,...)
%
% INPUTS:
%    model:             The COBRA  model structure. 
%    type:              the basic field to get annotations for (e.g. rxn, met, or
%                       model
%
% OPTIONAL INPUTS:
%    varagin:           Additional arguments as parameter/value pairs.
%                        * bioQualifiers - A Cell array of BioQualifiers to look for if not provided, or empty, all bioQualifiers defined in `getBioQualifiers()` will be used
%                        * ids - A specific ID or IDs to get the annotation data for. Cannot be combined with type model.(Default: model.([type 's'])). 
%
% OUTPUT:
%    annotations:       A struct array with the following structure:
%                        * .annotation -> A Struct array with one element per element of the given type.
%                        * .annotation.cvterms -> A Struct array of controlled vocabulary terms with one element per qualifier used
%                        * .annotation.cvterms.qualifier -> the bioQualifier for this all ressourcesof this cvterm 
%                        * .annotation.cvterms.qualifierType -> the Qualifier type (modelQualifier or bioQualifier) for all ressources of this cvterm 
%                        * .annotation.cvterms.ressources -> struct with the following fields:
%                        * .annotation.cvterms.ressources.database -> the database for this ressource.
%                        * .annotation.cvterms.ressources.id-> The ID in the database for this ressource.


[defaultBioQualifiers,standardQualifiers] = getBioQualifiers();
if ~strcmpi(type,'model')
    defaultIDs = model.([lower(type) 's']);
else
    defaultIDs = {'model'};
end

parser = inputParser();

parser.addParameter('bioQualifiers',defaultBioQualifiers,@(x) ischar(x) || iscell(x))
parser.addParameter('ids',defaultIDs,@(x) ischar(x) || iscell(x));

parser.parse(varargin{:});

bioQualifiers = parser.Results.bioQualifiers;
ids = parser.Results.ids;

% we have to handle some things special for model annotations (e.g. they can
% have model qualifiers)
if strcmp(type,'model')
    numElements = 1;
    annotations = cell(1);
    modelAnnot = true;
    bioQualifiers = [strcat('m',bioQualifiers),strcat('b',bioQualifiers)];    
    ids = {'model'};
else
    numElements = length(ids);
    modelAnnot = false;
    modelQualString = '';
    relPos = ismember(model.([type 's']),ids);
end

modelFields = fieldnames(model);
unusedFields = true(size(modelFields));
databaseFields = getDatabaseMappings(type);
% ok, these can be converted
databaseFields = databaseFields(ismember(databaseFields(:,3),modelFields),:);
% so we first convert the databaseFields
% but we will first filter so that we only get the standard qualifiers for
% each db.
databaseFields = databaseFields(ismember(databaseFields(:,2),standardQualifiers(:,1)),:);

% now we know the database fields, so we will check for the non default
% annotation fields
relfields = modelFields(cellfun(@(x) strncmp(x,type,length(type)),modelFields));        
annotationsFields = relfields(cellfun(@(x) any(cellfun(@(y) strncmp(x,[type, y],length([type y])),bioQualifiers)),relfields));

% now, we can initialize the result cell array. 
% this array has the following form:
% Dim1 -> one entry per element
% Dim2 -> one entry per either database or annotations Field
% Dim3 -> 1: qualifierType ; 2: qualifier ; 3: Database ; 4: IDs ; 
resultArray = cell(numElements,size(databaseFields,1)+ numel(annotationsFields), 4);
arrayDim2Pos = 1;
for i = 1:size(databaseFields,1)
    cField = databaseFields{i,3};
    cSource = databaseFields(i,1);
    cQual = databaseFields(i,2);
    if modelAnnot
        cValues = {model.(cField)};    
    else
        cValues = model.(cField)(relPos);
    end
    cQualType = databaseFields(i,6);
    resultArray(:,arrayDim2Pos,:) = [repmat(cQualType,numElements,1),repmat(cQual,numElements,1),repmat(cSource,numElements,1),cValues];
    arrayDim2Pos = arrayDim2Pos + 1;
end

for i = 1:numel(annotationsFields)
    cField = annotationsFields{i};
    if modelAnnot
        cValues = {model.(annotationsFields{i})};    
    else
        cValues = model.(annotationsFields{i})(relPos);
    end
    % get the correct qualifier for this field.
    for i = 1:numel(bioQualifiers)
        cQual = bioQualifiers{i};
        if strncmp(cField,[type, cQual],length([type cQual]))
            cField = cField((length([type cQual])+1):end);
            break;
        end
    end
    cQualType = 'bioQualifier';
    if modelAnnot
        if cQual(1) == 'm'
            cQualType = 'modelQualifier';
        end
        cQual = cQual(2:end);
    end
    cSource = convertSBMLID(regexprep(cField,'ID$',''),false);
    resultArray(:,arrayDim2Pos,:) = [repmat({cQualType},numElements,1),repmat({cQual},numElements,1),repmat({cSource},numElements,1),cValues];
    arrayDim2Pos = arrayDim2Pos + 1;
end

% now, build a struct out of this cell array.
ressourceStruct = struct('database','','id','');
ressourceStruct(1) = []; %Clear the empty element.
cvtermStruct = struct(struct('qualifier','','qualifierType','','ressources',ressourceStruct));
cvtermStruct(1) = [];
annotations = struct('id','rxn1','cvterms',cvtermStruct);
annotations(numElements).id = '';

% as a reminder, this is the structure of the array:
% Dim1 -> one entry per element
% Dim2 -> one entry per either database or annotations Field
% Dim3 -> 1: qualifierType ; 2: qualifier ; 3: Database ; 4: IDs ; 

for i = 1:numElements
    cvtermsIndex = 1;
    annotations(i).id = ids{i};
    currentCVtermsStruct = cvtermStruct;    
    % take all qualifierTypes which have non empty entries for this element
    relArray = resultArray(i,:,:);
    currentQualifiers = unique(resultArray(i,:,2));
    % we go over each qualifier, and will then look for the qualifierType to
    % distinguish
    for j = 1:numel(currentQualifiers)
        % get the current qualifier
        cQualifier = currentQualifiers{j};
        qualIndex = strcmp(resultArray(i,:,2),cQualifier);
        qualifierTypes = unique(resultArray(i,qualIndex,1));        
        for k = 1:numel(qualifierTypes)
            % and qualifier typee
            cQualifierType = qualifierTypes{k};            
            cStruct = struct('qualifier',cQualifier,'qualifierType',cQualifierType,'ressources',ressourceStruct);
            cDBArray = {};
            cIDArray = {};
            relIndices = strcmp(resultArray(i,:,1),cQualifierType) & qualIndex;
            databases = resultArray(i,relIndices,3);
            dbids = resultArray(i,relIndices,4);
            for cdb = 1:length(databases)
                % and the ressources
                cDatabase = databases{cdb};
                cIDs = strsplit(dbids{cdb},'; '); %IDs are split by ; in the fields.                
                if ~all(cellfun(@isempty, cIDs))
                    cDBArray = [cDBArray , repmat({cDatabase},1,numel(cIDs))];
                    cIDArray = [cIDArray , cIDs];
                end
            end
            if ~isempty(cIDArray)
                % if we have at least one ID
                cRessourcestruct = ressourceStruct;
                cRessourceStruct(numel(cIDArray)).id = '';
                [cRessourceStruct(:).id] = deal(cIDArray{:});
                [cRessourceStruct(:).database] = deal(cDBArray{:});
                cStruct.ressources = cRessourceStruct();           
                currentCVtermsStruct(cvtermsIndex) = cStruct;
                cvtermsIndex = cvtermsIndex + 1;                           
            end
            annotations(i).cvterms = currentCVtermsStruct; %Either empty or with data.
        end        
    end
end

end