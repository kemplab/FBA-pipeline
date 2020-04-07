
function [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_ccle(input_file,constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein)
    
    % load data - CCLE
    [~,~,raw] = xlsread(input_file,'Samples - CCLE');

    % samples - CCLE
    ccle_types = {};
    ccle_needs_radiation = 0;
    for row = 3:26
        value = raw(row,1);
        value = value{1};

        if ~isnan(value)
            if ischar(value)
                if strcmpi(value,'X')
                    switch row
                        case 3, ccle_types{end+1} = 'Autonomic ganglia';
                        case 4, ccle_types{end+1} = 'Biliary tract';
                        case 5, ccle_types{end+1} = 'Bone';
                        case 6, ccle_types{end+1} = 'Breast';
                        case 7, ccle_types{end+1} = 'Central nervous system';
                        case 8, ccle_types{end+1} = 'Endometrium';
                        case 9, ccle_types{end+1} = 'Oesophagus';
                        case 10, ccle_types{end+1} = 'Haematopoietic and lymphoid tissue';
                        case 11, ccle_types{end+1} = 'Kidney';
                        case 12, ccle_types{end+1} = 'Large intestine';
                        case 13, ccle_types{end+1} = 'Liver';
                        case 14, ccle_types{end+1} = 'Lung';
                        case 15, ccle_types{end+1} = 'Ovary';
                        case 16, ccle_types{end+1} = 'Pancreas';
                        case 17, ccle_types{end+1} = 'Pleura';
                        case 18, ccle_types{end+1} = 'Prostate';
                        case 19, ccle_types{end+1} = 'Salivary gland';
                        case 20, ccle_types{end+1} = 'Skin';
                        case 21, ccle_types{end+1} = 'Small intestine';
                        case 22, ccle_types{end+1} = 'Soft tissue';
                        case 23, ccle_types{end+1} = 'Stomach';
                        case 24, ccle_types{end+1} = 'Thyroid';
                        case 25, ccle_types{end+1} = 'Upper aerodigestive tract';
                        case 26, ccle_types{end+1} = 'Urinary tract';
                    end
                else
                    error('ERROR - Samples - CCLE - Values must either be X or blank')
                end
            else
                error('ERROR - Samples - CCLE - Values must either be X or blank')
            end
        end
    end
 
    % needs mutation data
    row = 3;
    value = raw(row,4);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                ccle_needs_mutation = 1;
            else
                error('ERROR - Samples - CCLE - Restrict samples to those with mutation data - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - CCLE - Restrict samples to those with mutation data - Value must either be X or blank')
        end
    elseif (use_mutation == 1) && (~isempty(samples_ccle))
        error('ERROR - Samples - CCLE - Restrict samples to those with mutation data - Must choose this option if using mutation data')
    else
        ccle_needs_mutation = 0;
    end

    % needs radiation response
    row = 4;
    value = raw(row,4);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                ccle_needs_radiation = 1;
            else
                error('ERROR - Samples - CCLE - Restrict samples to those with radiation response - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - CCLE - Restrict samples to those with radiation response - Value must either be X or blank')
        end
    else
        ccle_needs_radiation = 0;
    end

    % needs any drug response
    row = 5;
    value = raw(row,4);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                ccle_needs_any_drug = 1;
            else
                error('ERROR - Samples - CCLE - Restrict samples to those with any drug response - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - CCLE - Restrict samples to those with any drug response - Value must either be X or blank')
        end
    elseif (length(required_drugs) > 0) && (~isempty(samples_ccle))
        error('ERROR - Samples - CCLE - Restrict samples to those with mutation data - Must choose this option if using mutation data')
    else
        ccle_needs_any_drug = 0;
    end

    % get CCLE samples
    if ~isempty(ccle_types)

        % get list of all CCLE samples
        files = dir('../data/clinical/CCLE/');
        fn = {};
        for i = 3:length(files)
            temp = strsplit(files(i).name,'.csv');
            fn{end+1} = temp{1};
        end
        
        % get list of all CCLE vmax data
        files = dir('../data/vmax/CCLE/no_mutation/');
        fn_protein = {};
        for i = 3:length(files)
            temp = strsplit(files(i).name,'.csv');
            fn_protein{end+1} = temp{1};
        end
        
        % get all samples with mutation data
        files = dir('../data/mutation/CCLE/');
        fn_mutation = {};
        for i = 3:length(files)
            temp = strsplit(files(i).name,'.csv');
            fn_mutation{end+1} = temp{1};
        end
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
        fn_mutation = [fn_mutation,idh1_patients_H133Q,idh1_patients_A134D,idh1_patients_R100Q,idh1_patients_R132H,idh1_patients_R132C,idh1_patients_R132G,idh1_patients_R132W,idh1_patients_R132A,idh1_patients_R132Q,idh1_patients_R132K,idh1_patients_R132N];
        
        % iterate over samples
        for i = 1:length(fn)
            keep = 1;

            % load clinical file
            f = fopen(sprintf('../data/clinical/CCLE/%s.csv',fn{i}),'r');
            data = textscan(f,'%s %s','Delimiter',',');
            fclose(f);
            
            % determine if has vmax data
            if ~any(strcmp(fn_protein,fn{i}))
                keep = 0;

            % determine if in selected cancer type
            elseif ~any(strcmp(ccle_types,data{2}{find(strcmp(data{1},'PRIMARY SITE'))}))
                keep = 0;

            % if needs radiation response
            elseif strcmp(data{2}{find(strcmp(data{1},'RADIATION AUC'))},'') && (ccle_needs_radiation == 1)
                keep = 0;

            % if needs mutation data
            elseif (ccle_needs_mutation == 1) && ~any(strcmp(fn_mutation,fn{i}))
                keep = 0;

            % if needs any drug response 
            elseif ccle_needs_any_drug == 1     
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

            % add to samples list
            if keep == 1
                samples_all{end+1} = fn{i};
                samples_all_source{end+1} = 'CCLE';
                samples_all_source_protein{end+1} = 'CCLE';
            end
        end
    end
end