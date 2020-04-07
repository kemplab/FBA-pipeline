
function [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_tcga(input_file,constraints_proliferation,required_drugs,use_mutation)
    
    % load data - TCGA
    [~,~,raw] = xlsread(input_file,'Samples - TCGA');

    % samples - TCGA
    tcga_types = {};
    tcga_needs_proliferation = 0;
    tcga_needs_radiation = 0;
    for row = 3:35
        value = raw(row,1);
        value = value{1};

        if ~isnan(value)
            if ischar(value)
                if strcmpi(value,'X')
                    switch row
                        case 3, tcga_types{end+1} = 'ACC';
                        case 4, tcga_types{end+1} = 'BLCA';
                        case 5, tcga_types{end+1} = 'BRCA';
                        case 6, tcga_types{end+1} = 'CESC';
                        case 7, tcga_types{end+1} = 'CHOL';
                        case 8, tcga_types{end+1} = 'COAD';
                        case 9, tcga_types{end+1} = 'DLBC';
                        case 10, tcga_types{end+1} = 'ESCA';
                        case 11, tcga_types{end+1} = 'GBM';
                        case 12, tcga_types{end+1} = 'HNSC';
                        case 13, tcga_types{end+1} = 'KICH';
                        case 14, tcga_types{end+1} = 'KIRC';
                        case 15, tcga_types{end+1} = 'KIRP';
                        case 16, tcga_types{end+1} = 'LAML';
                        case 17, tcga_types{end+1} = 'LGG';
                        case 18, tcga_types{end+1} = 'LIHC';
                        case 19, tcga_types{end+1} = 'LUAD';
                        case 20, tcga_types{end+1} = 'LUSC';
                        case 21, tcga_types{end+1} = 'MESO';
                        case 22, tcga_types{end+1} = 'OV';
                        case 23, tcga_types{end+1} = 'PAAD';
                        case 24, tcga_types{end+1} = 'PCPG';
                        case 25, tcga_types{end+1} = 'PRAD';
                        case 26, tcga_types{end+1} = 'READ';
                        case 27, tcga_types{end+1} = 'SARC';
                        case 28, tcga_types{end+1} = 'SKCM';
                        case 29, tcga_types{end+1} = 'STAD';
                        case 30, tcga_types{end+1} = 'TGCT';
                        case 31, tcga_types{end+1} = 'THCA';
                        case 32, tcga_types{end+1} = 'THYM';
                        case 33, tcga_types{end+1} = 'UCEC';
                        case 34, tcga_types{end+1} = 'UCS';
                        case 35, tcga_types{end+1} = 'UVM';
                    end
                else
                    error('ERROR - Samples - TCGA - Values must either be X or blank')
                end
            else
                error('ERROR - Samples - TCGA - Values must either be X or blank')
            end
        end
    end
    
    % old or new protein abundance
    if length(tcga_types) > 0
        for row = 3:4
            value = raw(row,5);
            value = value{1};
            if ~isnan(value)
                if ischar(value)
                    if strcmpi(value,'X')
                        if exist('tcga_protein','var')
                            error('ERROR - Samples - TCGA - Old or new protein abundance - Can only choose one option')
                        else
                            if row == 3
                                tcga_protein = 'TCGA_old';
                            elseif row == 4
                                tcga_protein = 'TCGA';
                            end
                        end
                    else
                        error('ERROR - Samples - TCGA - Old or new protein abundance - Values must either be X or blank')
                    end
                else
                    error('ERROR - Samples - TCGA - Old or new protein abundance - Values must either be X or blank')
                end
            end
        end
        if ~exist('tcga_protein','var')
            error('ERROR - Samples - TCGA - Old or new protein abundance - Must select one option')
        end
    else
        tcga_protein = NaN;
    end
    
    % include primary tumor samples
    row = 8;
    value = raw(row,5);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                tcga_include_tumor = 1;
            else
                error('ERROR - Samples - TCGA - Include primary tumor samples - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - TCGA - Include primary tumor samples - Value must either be X or blank')
        end
    else
        tcga_include_tumor = 0;
    end

    % include normal tissue samples
    row = 9;
    value = raw(row,5);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                tcga_include_normal = 1;
            else
                error('ERROR - Samples - TCGA - Include normal tissue samples - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - TCGA - Include normal tissue samples - Value must either be X or blank')
        end
    else
        tcga_include_normal = 0;
    end

    % needs doubling time
    row = 13;
    value = raw(row,5);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                tcga_needs_proliferation = 1;
            else
                error('ERROR - Samples - TCGA - Restrict samples to those with doubling time - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - TCGA - Restrict samples to those with doubling time - Value must either be X or blank')
        end
    elseif (constraints_proliferation == 1) && (~isempty(samples_tcga))
        error('ERROR - Samples - TCGA - Restrict samples to those with doubling time - Must choose this option if imposing proliferation constraint')
    else
        tcga_needs_proliferation = 0;
    end

    % needs mutation data
    row = 14;
    value = raw(row,5);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                tcga_needs_mutation = 1;
            else
                error('ERROR - Samples - TCGA - Restrict samples to those with mutation data - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - TCGA - Restrict samples to those with mutation data - Value must either be X or blank')
        end
    elseif (use_mutation == 1) && (~isempty(samples_tcga))
        error('ERROR - Samples - TCGA - Restrict samples to those with mutation data - Must choose this option if using mutation data')
    else
        tcga_needs_mutation = 0;
    end

    % needs radiation response
    row = 15;
    value = raw(row,5);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                tcga_needs_radiation = 1;
            else
                error('ERROR - Samples - TCGA - Restrict samples to those with radiation response - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - TCGA - Restrict samples to those with radiation response - Value must either be X or blank')
        end
    else
        tcga_needs_radiation = 0;
    end

    % needs any drug response
    row = 16;
    value = raw(row,5);
    value = value{1};

    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                tcga_needs_any_drug = 1;
            else
                error('ERROR - Samples - TCGA - Restrict samples to those with any drug response - Value must either be X or blank')
            end
        else
            error('ERROR - Samples - TCGA - Restrict samples to those with any drug response - Value must either be X or blank')
        end
    elseif (length(required_drugs) > 0) && (~isempty(samples_tcga))
        error('ERROR - Samples - TCGA - Restrict samples to those with mutation data - Must choose this option if using mutation data')
    else
        tcga_needs_any_drug = 0;
    end

    % get TCGA samples
    samples_all = {};
    samples_all_source = {};
    samples_all_source_protein = {};
    if ~isempty(tcga_types)

        % get list of all TCGA samples
        files = dir('../data/clinical/TCGA/');
        fn = {};
        for i = 3:length(files)
            temp = strsplit(files(i).name,'.csv');
            fn{end+1} = temp{1};
        end
        
        % get list of all TCGA vmax data
        files = dir(sprintf('../data/vmax/%s/no_mutation/',tcga_protein));
        fn_protein = {};
        for i = 3:length(files)
            temp = strsplit(files(i).name,'.csv');
            fn_protein{end+1} = temp{1};
        end
        
        % get all samples with mutation data
        files = dir('../data/mutation/TCGA/');
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
            f = fopen(sprintf('../data/clinical/TCGA/%s.csv',fn{i}),'r');
            data = textscan(f,'%s %s','Delimiter',',');
            fclose(f);
            
            % determine if has vmax data
            if ~any(strcmp(fn_protein,fn{i}))
                keep = 0;

            % determine if in selected cancer type
            elseif ~any(strcmp(tcga_types,data{2}{find(strcmp(data{1},'COHORT'))}))
                keep = 0;

            % if primary tumor
            elseif strcmp(data{2}{find(strcmp(data{1},'TYPE'))},'TUMOR') && (tcga_include_tumor == 0)
                keep = 0;

            % if normal tissue
            elseif strcmp(data{2}{find(strcmp(data{1},'TYPE'))},'NORMAL') && (tcga_include_normal == 0)
                keep = 0;

            % if needs doubling time
            else
                [~, status] = str2num(data{2}{find(strcmp(data{1},'PROLIFERATION [1/hr]'))});
                if (status == 0) && (tcga_needs_proliferation == 1)
                    keep = 0;

                % if needs radiation response
                elseif strcmp(data{2}{find(strcmp(data{1},'RESPONSE RADIATION'))},'') && (tcga_needs_radiation == 1)
                    keep = 0;

                % if needs mutation data
                elseif (tcga_needs_mutation == 1) && ~any(strcmp(fn_mutation,fn{i}))
                    keep = 0;

                % if needs any drug response 
                elseif tcga_needs_any_drug == 1     
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
                samples_all_source{end+1} = 'TCGA';
                samples_all_source_protein{end+1} = tcga_protein;
            end
        end
    end
end