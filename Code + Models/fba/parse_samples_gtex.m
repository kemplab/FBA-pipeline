
function [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_gtex(input_file,constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein)
    
    % load data - GTEx
    [~,~,raw] = xlsread(input_file,'Samples - GTEx');

    % samples - GTEx
    gtex_types = {};
    for row = 3:56
        value = raw(row,1);
        value = value{1};

        if ~isnan(value)
            if ischar(value)
                if strcmpi(value,'X')
                    switch row
                        case 3, gtex_types{end+1} = 'Adipose - Subcutaneous';
                        case 4, gtex_types{end+1} = 'Adipose - Subcutaneous';
                        case 5, gtex_types{end+1} = 'Adrenal Gland';
                        case 6, gtex_types{end+1} = 'Artery - Aorta';
                        case 7, gtex_types{end+1} = 'Artery - Coronary';
                        case 8, gtex_types{end+1} = 'Artery - Tibial';
                        case 9, gtex_types{end+1} = 'Bladder';
                        case 10, gtex_types{end+1} = 'Brain - Amygdala';
                        case 11, gtex_types{end+1} = 'Brain - Anterior cingulate cortex (BA24)';
                        case 12, gtex_types{end+1} = 'Brain - Caudate (basal ganglia)';
                        case 13, gtex_types{end+1} = 'Brain - Cerebellar Hemisphere';
                        case 14, gtex_types{end+1} = 'Brain - Cerebellum';
                        case 15, gtex_types{end+1} = 'Brain - Cortex';
                        case 16, gtex_types{end+1} = 'Brain - Frontal Cortex (BA9)';
                        case 17, gtex_types{end+1} = 'Brain - Hippocampus';
                        case 18, gtex_types{end+1} = 'Brain - Hypothalamus';
                        case 19, gtex_types{end+1} = 'Brain - Nucleus accumbens (basal ganglia)';
                        case 20, gtex_types{end+1} = 'Brain - Putamen (basal ganglia)';
                        case 21, gtex_types{end+1} = 'Brain - Spinal cord (cervical c-1)';
                        case 22, gtex_types{end+1} = 'Brain - Substantia nigra';
                        case 23, gtex_types{end+1} = 'Breast - Mammary Tissue';
                        case 24, gtex_types{end+1} = 'Cells - EBV-transformed lymphocytes';
                        case 25, gtex_types{end+1} = 'Cells - Leukemia cell line (CML)';
                        case 26, gtex_types{end+1} = 'Cells - Transformed fibroblasts';
                        case 27, gtex_types{end+1} = 'Cervix - Ectocervix';
                        case 28, gtex_types{end+1} = 'Cervix - Endocervix';
                        case 29, gtex_types{end+1} = 'Colon - Sigmoid';
                        case 30, gtex_types{end+1} = 'Colon - Transverse';
                        case 31, gtex_types{end+1} = 'Esophagus - Gastroesophageal Junction';
                        case 32, gtex_types{end+1} = 'Esophagus - Mucosa';
                        case 33, gtex_types{end+1} = 'Esophagus - Muscularis';
                        case 34, gtex_types{end+1} = 'Fallopian Tube';
                        case 35, gtex_types{end+1} = 'Heart - Atrial Appendage';
                        case 36, gtex_types{end+1} = 'Heart - Left Ventricle';
                        case 37, gtex_types{end+1} = 'Kidney - Cortex';
                        case 38, gtex_types{end+1} = 'Liver';
                        case 39, gtex_types{end+1} = 'Lung';
                        case 40, gtex_types{end+1} = 'Minor Salivary Gland';
                        case 41, gtex_types{end+1} = 'Muscle - Skeletal';
                        case 42, gtex_types{end+1} = 'Nerve - Tibial';
                        case 43, gtex_types{end+1} = 'Ovary';
                        case 44, gtex_types{end+1} = 'Pancreas';
                        case 45, gtex_types{end+1} = 'Pituitary';
                        case 46, gtex_types{end+1} = 'Prostate';
                        case 47, gtex_types{end+1} = 'Skin - Not Sun Exposed (Suprapubic)';
                        case 48, gtex_types{end+1} = 'Skin - Sun Exposed (Lower leg)';
                        case 49, gtex_types{end+1} = 'Small Intestine - Terminal Ileum';
                        case 50, gtex_types{end+1} = 'Spleen';
                        case 51, gtex_types{end+1} = 'Stomach';
                        case 52, gtex_types{end+1} = 'Testis';
                        case 53, gtex_types{end+1} = 'Thyroid';
                        case 54, gtex_types{end+1} = 'Uterus';
                        case 55, gtex_types{end+1} = 'Vagina';
                        case 56, gtex_types{end+1} = 'Whole Blood';
                    end
                else
                    error('ERROR - Samples - GTEx - Values must either be X or blank')
                end
            else
                error('ERROR - Samples - GTEx - Values must either be X or blank')
            end
        end
    end
 
    % get GTEx samples
    if ~isempty(gtex_types)

        % get list of all GTEx samples
        files = dir('../data/clinical/GTEx/');
        fn = {};
        for i = 3:length(files)
            temp = strsplit(files(i).name,'.csv');
            fn{end+1} = temp{1};
        end
        
        % get list of all GTEx vmax data
        files = dir('../data/vmax/GTEx/no_mutation/');
        fn_protein = {};
        for i = 3:length(files)
            temp = strsplit(files(i).name,'.csv');
            fn_protein{end+1} = temp{1};
        end
        
        % iterate over samples
        for i = 1:length(fn)
            keep = 1;

            % load clinical file
            f = fopen(sprintf('../data/clinical/GTEx/%s.csv',fn{i}),'r');
            data = textscan(f,'%s %s','Delimiter',',');
            fclose(f);
            
            % determine if has vmax data
            if ~any(strcmp(fn_protein,fn{i}))
                keep = 0;

            % determine if in selected cancer type
            elseif ~any(strcmp(gtex_types,data{2}{find(strcmp(data{1},'TISSUE TYPE'))}))
                keep = 0;

            % required drugs
            elseif length(required_drugs) > 0
                keep = 0;
            end

            % add to samples list
            if keep == 1
                samples_all{end+1} = fn{i};
                samples_all_source{end+1} = 'GTEx';
                samples_all_source_protein{end+1} = 'GTEx';
            end
        end
    end
end