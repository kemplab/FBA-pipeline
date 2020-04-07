
function [vmax_forward,vmax_reverse] = calculate_vmax(model,samples_all,samples_all_source,samples_all_source_protein,use_mutation,knockdown_genes,knockdown_fractions,gene_knockdown_screen)

    % initialize vmax
    vmax_forward = nan(length(model.rxns),length(samples_all));
    vmax_reverse = nan(length(model.rxns),length(samples_all));
    
    % if gene knockdown, parse GPR rules and reaction kcat's
    if (length(knockdown_genes) > 0)
 
        % parse GPR rules
        parsedGPR = GPRparser(model);
                    
        % load forward reaction kcat's
        f = fopen('../data/kcat/kcat_forward.csv','r');
        kcat_forward_original = textscan(f,'%s %f %s','Delimiter',',','headerLines',1);
        fclose(f);

        % load reverse reaction kcat's
        f = fopen('../data/kcat/kcat_reverse.csv','r');
        kcat_reverse_original = textscan(f,'%s %f %s','Delimiter',',','headerLines',1);
        fclose(f);
    end
    
    % load list of TCGA patients with IDH1 mutations
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
    idh1_patients_all = [idh1_patients_H133Q,idh1_patients_A134D,idh1_patients_R100Q,idh1_patients_R132H,idh1_patients_R132C,idh1_patients_R132G,idh1_patients_R132W,idh1_patients_R132A,idh1_patients_R132Q,idh1_patients_R132K,idh1_patients_R132N];

    % IDH1 kcat, forward reaction
    idh1_normal_H133Q = 45*60*60;
    idh1_normal_R132A = 10.4*60*60;
    idh1_normal_R132G = 9.3*60*60;
    idh1_normal_R132Q = 9.2*60*60;
    idh1_normal_R132K = 7.2*60*60;
    idh1_normal_R100Q = 5.6*60*60;
    idh1_normal_R132C = 4.4*60*60;
    idh1_normal_R132H = 2.4*60*60;
    idh1_normal_A134D = 2.3*60*60;
    idh1_normal_R132W = 1.21*60*60;
    idh1_normal_R132N = 0.047*60*60;

    % IDH1 kcat, neomorphic reaction
    idh1_neomorphic_H133Q = 0*60*60;
    idh1_neomorphic_A134D = 0*60*60;
    idh1_neomorphic_R100Q = 0.34*60*60;
    idh1_neomorphic_R132A = 0.37*60*60;
    idh1_neomorphic_R132W = 0.54*60*60;
    idh1_neomorphic_R132K = 0.57*60*60;
    idh1_neomorphic_R132N = 0.79*60*60;
    idh1_neomorphic_R132G = 1.59*60*60;
    idh1_neomorphic_R132C = 1.60*60*60;
    idh1_neomorphic_R132H = 4.2*60*60;
    idh1_neomorphic_R132Q = 4.7*60*60;

    % iterate over samples
    for i = 1:length(samples_all)

        % load vmax data
        if (use_mutation == 1) || ( ( (exist(sprintf('../data/mutation/%s/%s.csv',samples_all_source{i},samples_all{i}),'file')) || (any(strcmp(idh1_patients_all,samples_all{i}))) ) && (use_mutation == 2))
            f = fopen(sprintf('../data/vmax/%s/yes_mutation/%s.csv',samples_all_source_protein{i},samples_all{i}),'r');
        else
            f = fopen(sprintf('../data/vmax/%s/no_mutation/%s.csv',samples_all_source_protein{i},samples_all{i}),'r');
        end
        vmax = textscan(f,'%s %f %f','Delimiter',',','HeaderLines',1);
        vmax{1} = strrep(vmax{1},'"','');
        fclose(f);
        
        % put in vmax values
        vmax_forward(1:length(vmax{1}),i) = vmax{2};
        vmax_reverse(1:length(vmax{1}),i) = vmax{3};
 
        % if gene knockdown, change vmax values for reactions containing that gene
        if (length(knockdown_genes) > 0)
             
            % load protein expression
            f = fopen(sprintf('../data/protein/output/%s/%s.csv',samples_all_source_protein{i},samples_all{i}),'r');
            protein = textscan(f,'%s %f','Delimiter',',');
            protein{1} = strrep(protein{1},'"','');
            fclose(f);
            
            % if using mutation data
            if ( (exist(sprintf('../data/mutation/%s/%s.csv',samples_all_source{i},samples_all{i}),'file')) && (use_mutation == 1) ) || ( (exist(sprintf('../data/mutation/%s/%s.csv',samples_all_source{i},samples_all{i}),'file')) && (use_mutation == 2) )

                % load mutation data
                f = fopen(sprintf('../data/mutation/%s/%s.csv',samples_all_source{i},samples_all{i}),'r');
                mutation = textscan(f,'%s %f','Delimiter',',','headerLines',1);
                fclose(f);

                % multiply protein expression data by envision scores
                for j = 1:length(mutation{1})
                    if any(strcmp(protein{1},mutation{1}{j}))
                        protein{2}(find(strcmp(protein{1},mutation{1}{j}))) = protein{2}(find(strcmp(protein{1},mutation{1}{j}))) * mutation{2}(j);
                    end
                end
            end
            
            % gene knockdowns
            for j = 1:length(knockdown_genes)
                if any(strcmp(protein{1},model.geneSymbols{find(strcmp(model.genes,knockdown_genes{j}))}))
                    protein{2}(find(strcmp(protein{1},model.geneSymbols{find(strcmp(model.genes,knockdown_genes{j}))}))) = protein{2}(find(strcmp(protein{1},model.geneSymbols{find(strcmp(model.genes,knockdown_genes{j}))}))) * knockdown_fractions(j);
                end
            end
            
            % convert ppm to mmol/gDW
            m = 50; % average mass of protein, in kDa
            fdw = 0.5; % fraction of dry weight that is protein
            Na = 6.0221409*10^23; % Avogadro constant
            gDa = 1.66054*10^-24; % number of grams per dalton

            protein{2} = protein{2}/1000000; % protein/million --> protein/protein
            protein{2} = protein{2} / Na; % protein/protein --> mol/protein
            protein{2} = protein{2} * 1000; % mol/protein --> mmol/protein
            protein{2} = protein{2} / (m*1000); % mmol/protein --> mmol/Daltons_protein
            protein{2} = protein{2} / gDa; % mmol/Daltons --> mmol/grams_protein
            protein{2} = protein{2} * fdw; % mmol/grams_protein --> mmol/gDW
            
            % get reactions that gene(s) is(are) involved in
            reaction_list = [];
            for j = 1:length(knockdown_genes)
                reaction_list = [reaction_list,find(full(model.rxnGeneMat(:,strcmp(model.genes,knockdown_genes{j}))))'];
            end
            reaction_list = unique(reaction_list);
            
            % get protein expression for each reaction with grRule and EC number containing gene(s) being knockdown
            protein_rxn = nan(1,length(model.rxns));
            for j = reaction_list
                if (~isempty(model.rules{j})) && (~isempty(model.rxnECNumbers{j}))

                    % get genes associated with reaction
                    genes = find(model.rxnGeneMat(j,:));
                    allfound = 1;
                    for k = 1:length(genes)
                        if ~any(strcmp(protein{1},model.geneSymbols{genes(k)}));
                            allfound = 0;
                            break
                        end
                    end
                    if allfound == 1

                        % create model structure for selectGeneFromGPR
                        modelRxn = {};
                        modelRxn.rxns = {model.rxns{j}};
                        modelRxn.genes = {};
                        for k = 1:length(genes)
                            modelRxn.genes{end+1} = model.genes{genes(k)};
                        end
                        modelRxn.genes = modelRxn.genes';

                        % calculate reaction expression
                        values = [];
                        for k = 1:length(genes)
                            values(end+1) = protein{2}(strcmp(protein{1},model.geneSymbols{genes(k)}));
                        end
                        % updated version of selectGeneFromGPR
                        protein_rxn(j) = selectGeneFromGPR(modelRxn,modelRxn.genes,values',{parsedGPR{j}});

                    else
                       protein_rxn(j) = NaN; 
                    end

                else
                    protein_rxn(j) = NaN;
                end
            end
            
            % load kcat values
            kcat_forward = kcat_forward_original;
            kcat_reverse = kcat_reverse_original;
            
            % alter kcat reactions for IDH1 normal reactions
            if (use_mutation == 1) or (use_mutation == 2)
                if any(strcmp(idh1_patients_R132N,samples_all{i}))
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132N;
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132N;
                    kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132N;
                elseif any(strcmp(idh1_patients_R132W,samples_all{i}))
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132W;
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132W;
                    kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132W;
                elseif any(strcmp(idh1_patients_A134D,samples_all{i}))
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_A134D;
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_A134D;
                    kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_A134D;
                elseif any(strcmp(idh1_patients_R132H,samples_all{i}))
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132H;
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132H;
                    kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132H;
                elseif any(strcmp(idh1_patients_R132C,samples_all{i}))
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132C;
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132C;
                    kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132C;
                elseif any(strcmp(idh1_patients_R100Q,samples_all{i}))
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R100Q;
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R100Q;
                    kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R100Q;
                elseif any(strcmp(idh1_patients_R132K,samples_all{i}))
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132K;
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132K;
                    kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132K;
                elseif any(strcmp(idh1_patients_R132Q,samples_all{i}))
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132Q;
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132Q;
                    kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132Q;
                elseif any(strcmp(idh1_patients_R132G,samples_all{i}))
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132G;
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132G;
                    kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132G;
                elseif any(strcmp(idh1_patients_R132A,samples_all{i}))
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132A;
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132A;
                    kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132A;
                elseif any(strcmp(idh1_patients_H133Q,samples_all{i}))
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_H133Q;
                    kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_H133Q;
                    kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_H133Q;
                end
            end
            
            % change vmax for reactions containing gene(s) being knockdown
            for j = reaction_list
                if ~strcmp(model.subSystems{j},'Custom Constraint') && ~strcmp(model.subSystems{j},'Custom Objective');

                    % forward vmax
                    if model.ub(j) > 0
                        if (~(isnan(protein_rxn(j)))) && (~(isnan(kcat_forward{2}(strcmp(kcat_forward{1},model.rxns{j})))))
                            vmax_forward(j,i) = protein_rxn(j) * kcat_forward{2}(strcmp(kcat_forward{1},model.rxns{j}));
                        end
                    end

                    % reverse vmax
                    if model.lb(j) < 0
                        if (~(isnan(protein_rxn(j)))) && (~(isnan(kcat_reverse{2}(strcmp(kcat_reverse{1},model.rxns{j})))))
                            vmax_reverse(j,i) = protein_rxn(j) * kcat_reverse{2}(strcmp(kcat_reverse{1},model.rxns{j}));
                        end
                    end
                end
            end
        end
    end
end