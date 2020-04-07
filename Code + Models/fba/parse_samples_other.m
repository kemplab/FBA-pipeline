
function [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_other(input_file,constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein)
    
    % load data - Other samples
    [~,~,raw] = xlsread(input_file,'Samples - Other');
    
    % other samples
    other_samples_folder = {};
    other_samples = {};

    for row = 2:size(raw,1)
        added = false;

        % folder name
        value = raw(row,1);
        value = value{1};
        if ~isnan(value)
            if ischar(value)
                
                % sample folder
                other_samples_folder{end+1} = value;
                
                % sample name
                value = raw(row,2);
                value = value{1};
                if ~isnan(value)
                    if ischar(value)
                        other_samples{end+1} = value;
                    elseif ischar(num2str(value))
                        other_samples{end+1} = num2str(value);
                    else
                        error('ERROR - Samples - Other - Sample name must be a valid string')
                    end
                else
                    error('ERROR - Samples - Other - Must enter sample name if provide sample folder')
                end
            end
        end
    end
    
    % needs doubling time
    row = 6;
    value = raw(row,4);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                other_needs_proliferation = 1;
            else
                error('ERROR - Samples - Other - Restrict samples to those with doubling time - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - Other - Restrict samples to those with doubling time - Value must either be X or blank')
        end
    elseif (constraints_proliferation == 1) && (~isempty(other_samples))
        error('ERROR - Samples - Other - Restrict samples to those with doubling time - Must choose this option if imposing proliferation constraint')
    else
        other_needs_proliferation = 0;
    end

    % needs mutation data
    row = 7;
    value = raw(row,4);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                other_needs_mutation = 1;
            else
                error('ERROR - Samples - Other - Restrict samples to those with mutation data - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - Other - Restrict samples to those with mutation data - Value must either be X or blank')
        end
    elseif (use_mutation == 1) && (~isempty(other_samples))
        error('ERROR - Samples - Other - Restrict samples to those with mutation data - Must choose this option if using mutation data')
    else
        other_needs_mutation = 0;
    end

    % needs any drug response
    row = 8;
    value = raw(row,4);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                other_needs_any_drug = 1;
            else
                error('ERROR - Samples - Other - Restrict samples to those with any drug response - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - Other - Restrict samples to those with any drug response - Value must either be X or blank')
        end
    elseif (length(required_drugs) > 0) && (~isempty(other_samples))
        error('ERROR - Samples - Other - Restrict samples to those with mutation data - Must choose this option if using mutation data')
    else
        other_needs_any_drug = 0;
    end
    
    % extract samples from _ALL_ folders
    remove = [];
    for i = 1:length(other_samples_folder)
        if strcmp(other_samples{i},'_ALL_')
            
            % remove index
            remove(end+1) = i;
            
            % get samples in folder
            files = dir(sprintf('../data/vmax/%s/no_mutation/',other_samples_folder{i}));
            fn_protein = {};
            for j = 3:length(files)
                temp = strsplit(files(j).name,'.csv');
                fn_protein{end+1} = temp{1};
            end
            
            % add samples
            for j = 1:length(fn_protein)
                other_samples_folder{end+1} = other_samples_folder{i};
                other_samples{end+1} = fn_protein{j};
            end 
        end
    end
    for i = fliplr(sort(remove))
        other_samples_folder(i) = [];
        other_samples(i) = [];
    end
    
    % patients with idh1 mutations
    idh1_patients_H133Q = importdata('../data/idh1/H133Q.txt')';
    idh1_patients_A134D = importdata('../data/idh1/A134D.txt')';
    idh1_patients_R100Q = importdata('../data/idh1/R100Q.txt')';
    idh1_patients_R132H = importdata('../data/idh1/R132H.txt')';
    idh1_patients_R132C = importdata('../data/idh1/R132C.txt')';
    idh1_patients_R132G = importdata('../data/idh1/R132G.txt')';
    idh1_patients_R132W = importdata('../data/idh1/R132W.txt')';
    idh1_patients_R132A = importdata('../data/idh1/R132A.txt')';
    idh1_patients_R132Q = importdata('../data/idh1/R132Q.txt')';
    idh1_patients_R132K = importdata('../data/idh1/R132K.txt')';
    idh1_patients_R132N = importdata('../data/idh1/R132N.txt')';
    
    % iterate over sample folders
    for i = 1:length(other_samples_folder)
        keep = 1;
        
        % folder
        temp = strsplit(other_samples_folder{i},'_');
        prefix = temp{1};
        switch prefix
            case 'TCGA'
                nonprotein_folder = 'TCGA';
            case 'NCI60'
                nonprotein_folder = 'NCI60';
            otherwise
                nonprotein_folder = other_samples_folder{i};
        end

        % check vmax data
        if ~(exist(sprintf('../data/vmax/%s',other_samples_folder{i})) == 7)
            error('ERROR - Samples - Other - Folder "%s" does not exist in vmax data',other_samples_folder{i})
        else
            files = dir(sprintf('../data/vmax/%s/no_mutation/',other_samples_folder{i}));
            fn_protein = {};
            for j = 3:length(files)
                temp = strsplit(files(j).name,'.csv');
                fn_protein{end+1} = temp{1};
            end
            if ~any(strcmp(fn_protein,other_samples{i}))
                error('ERROR - Samples - Other - Sample "%s" does not exist in folder "%s" of vmax data',other_samples{i},other_samples_folder{i})
            else
                
                % check mutation data
                if (other_needs_mutation == 1)
                    if ~(exist(sprintf('../data/mutation/%s',nonprotein_folder)) == 7)
                        keep = 0;
                    else
                        files = dir(sprintf('../data/mutation/%s/',nonprotein_folder));
                        fn_mutation = {};
                        for j = 3:length(files)
                            temp = strsplit(files(j).name,'.csv');
                            fn_mutation{end+1} = temp{1};
                        end
                        if strcmp(prefix,'TCGA')
                            fn_mutation = [fn_mutation,idh1_patients_H133Q,idh1_patients_A134D,idh1_patients_R100Q,idh1_patients_R132H,idh1_patients_R132C,idh1_patients_R132G,idh1_patients_R132W,idh1_patients_R132A,idh1_patients_R132Q,idh1_patients_R132K,idh1_patients_R132N];
                        end
                        if ~any(strcmp(fn_mutation,other_samples{i}))
                            keep = 0;
                        else
                
                            % check clinical data
                            if (other_needs_proliferation == 1) || (other_nees_any_drug == 1)
                                if ~(exist(sprintf('../data/clinical/%s',nonprotein_folder)) == 7)
                                    keep = 0;
                                else
                                    files = dir(sprintf('../data/clinical/%s/',nonprotein_folder));
                                    fn_protein = {};
                                    for j = 3:length(files)
                                        temp = strsplit(files(j).name,'.csv');
                                        fn_protein{end+1} = temp{1};
                                    end
                                    if ~any(strcmp(fn_protein,other_samples{i}))
                                        keep = 0;
                                    else
                                    
                                        % load clinical data
                                        f = fopen(sprintf('../data/clinical/%s/%s.csv',nonprotein_folder,other_samples{i}),'r');
                                        data = textscan(f,'%s %s','Delimiter',',');
                                        fclose(f);

                                        % if needs doubling time
                                        [~, status] = str2num(data{2}{find(strcmp(data{1},'PROLIFERATION [1/hr]'))});
                                        if (status == 0) && (other_needs_proliferation == 1)
                                            keep = 0;

                                        % if needs any drug response 
                                        elseif other_needs_any_drug == 1     
                                            response = 0;
                                            for j = 1:length(data{1})
                                                temp = strsplit(data{1}{j},' ');
                                                if length(temp) > 1
                                                    if strcmp(temp{2},'DRUG') && (~strcmp(data{2}{j},''))
                                                        response = 1;
                                                    end
                                                end
                                            end
                                            if response == 0
                                                keep = 0;
                                            end
                                            
                                            % required drugs
                                            if keep == 1
                                                for j = 1:length(required_drugs)

                                                    % is drug response available
                                                    if ~any(strcmp(data{1},sprintf('RESPONSE DRUG %s',required_drugs{j})))
                                                        keep = 0;
                                                    elseif strcmp(data{2}{find(strcmp(data{1},sprintf('RESPONSE DRUG %s',required_drugs{j})))},'')
                                                        keep = 0;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        % add to samples list
        if keep == 1
            samples_all{end+1} = other_samples{i};
            samples_all_source{end+1} = nonprotein_folder;
            samples_all_source_protein{end+1} = other_samples_folder{i};
        end
    end
end