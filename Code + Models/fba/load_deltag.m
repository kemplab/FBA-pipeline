
function [deltag_lower,deltag_upper] = load_deltag(model,recalculate_thermodynamics,default_concentration_lb,default_concentration_ub,concentration_ranges_all,concentration_ranges_sample_folder,concentration_ranges_sample,concentration_mets,concentration_values,samples_all,samples_all_source)

    % if recalculating thermodynamics
    if recalculate_thermodynamics == 1
    
        % load datasets for all samples
        conc_all_met = {};
        conc_all_lower = [];
        conc_all_upper = [];
        for i = 1:length(concentration_ranges_all)
        
            % load data
            f = fopen(sprintf('../data/concentration/_ALL_/%s.csv',concentration_ranges_all{i}),'r');
            data = textscan(f,'%s %f %f','Delimiter',',','headerLines',1);
            fclose(f);
            
            % add metabolite id, lower and upper bound
            for j = 1:length(data{1})
                
                % if metabolite already added
                if any(strcmp(conc_all_met,data{1}{j}))
                
                    % change lower bound if lower than original
                    if data{2}(j) < conc_all_lower(strcmp(conc_all_met,data{1}{j}))
                        conc_all_lower(strcmp(conc_all_met,data{1}{j})) = data{2}(j);
                    end
                    
                    % change upper bound if greater than original
                    if data{3}(j) > conc_all_upper(strcmp(conc_all_met,data{1}{j}))
                        conc_all_upper(strcmp(conc_all_met,data{1}{j})) = data{3}(j);
                    end  
                
                % if metabolite not already added
                else
                
                    % add info
                    conc_all_met{end+1} = data{1}{j};
                    conc_all_lower(end+1) = data{2}(j);
                    conc_all_upper(end+1) = data{3}(j);        
                end
            end
        end
        
        % load datasets for individual samples
        conc_individual_sample = {};
        conc_individual_dataset = {};
        conc_individual_met = {};
        conc_individual_lower = {};
        conc_individual_upper = {};
        for i = 1:length(concentration_ranges_sample)
        
            % add "all" concentrations to individual sample list
            conc_individual_sample{end+1} = concentration_ranges_sample{i};
            conc_individual_dataset{end+1} = concentration_ranges_sample_folder{i};
            conc_individual_met{end+1} = conc_all_met;
            conc_individual_lower{end+1} = conc_all_lower;
            conc_individual_upper{end+1} = conc_all_upper;
                
            % load data
            f = fopen(sprintf('../data/concentration/%s/%s.csv',concentration_ranges_sample_folder{i},concentration_ranges_sample{i}),'r');
            data = textscan(f,'%s %f %f','Delimiter',',','headerLines',1);
            fclose(f);
                
            % add metabolite id, lower and upper bound
            for j = 1:length(data{1})

                % if metabolite already added
                if any(strcmp(conc_individual_met{end},data{1}{j}))

                    % change lower bound if lower than original
                    if data{2}(j) < conc_individual_lower{end}(strcmp(conc_individual_met{end},data{1}{j}))
                        conc_individual_lower{end}(strcmp(conc_individual_met{end},data{1}{j})) = data{2}(j);
                    end

                    % change upper bound if greater than original
                    if data{3}(j) > conc_individual_upper{end}(strcmp(conc_individual_met{end},data{1}{j}))
                        conc_individual_upper{end}(strcmp(conc_individual_met{end},data{1}{j})) = data{3}(j);
                    end    

                % if metabolite not already added
                else

                    % add info
                    conc_individual_met{end}{end+1} = data{1}{j};
                    conc_individual_lower{end}(end+1) = data{2}(j);
                    conc_individual_upper{end}(end+1) = data{3}(j);
                end
            end
        end
        
        % load deltaG file
        f = fopen('../data/deltag/deltag.csv','r');
        data = textscan(f,'%s %f %f %f %f %f','Delimiter',',','headerLines',1);
        fclose(f);
        
        % calculate deltag ranges based on concentrations
        deltag_lower = -Inf*ones(length(model.rxns),length(samples_all));
        deltag_upper = Inf*ones(length(model.rxns),length(samples_all));
        
        % implement concentrations for all samples
        for i = 1:length(model.rxns)
            
            % if deltag data available
            if any(strcmp(data{1},model.rxns{i}))
            
                % get concentrations of reactants and products
                reactant_power = [];
                reactant_lower = [];
                reactant_upper = [];
                product_power = [];
                product_lower = [];
                product_upper = [];

                % iterate over reactants
                for j = find(model.S(:,i) < 0)'
                    reactant_power = full(model.S(j,i));

                    % if metabolite found in concentration data
                    if any(strcmp(conc_all_met,model.mets{j}))
                        reactant_lower(end+1) = conc_all_lower(strcmp(conc_all_met,model.mets{j}));
                        reactant_upper(end+1) = conc_all_upper(strcmp(conc_all_met,model.mets{j}));

                    % if metabolite not found in concentration data
                    else
                        reactant_lower(end+1) = default_concentration_lb;
                        reactant_upper(end+1) = default_concentration_ub;
                    end
                end

                % iterate over products
                for j = find(model.S(:,i) > 0)'
                    product_power = full(model.S(j,i));

                    % if metabolite found in concentration data
                    if any(strcmp(conc_all_met,model.mets{j}))
                        reactant_lower(end+1) = conc_all_lower(strcmp(conc_all_met,model.mets{j}));
                        reactant_upper(end+1) = conc_all_upper(strcmp(conc_all_met,model.mets{j}));

                    % if metabolite not found in concentration data
                    else
                        reactant_lower(end+1) = default_concentration_lb;
                        reactant_upper(end+1) = default_concentration_ub;
                    end
                end

                % calculate deltag bounds
                Q_max = prod(product_upper.^product_power)/prod(reactant_lower.^reactant_power);
                Q_min = prod(product_lower.^product_power)/prod(reactant_upper.^reactant_power);
                deltag_lower(i,:) = data{4}(strcmp(data{1},model.rxns{i})) - 8.3144598/1000*310.15*log(Q_max);
                deltag_upper(i,:) = data{4}(strcmp(data{1},model.rxns{i})) - 8.3144598/1000*310.15*log(Q_min);
            end
        end
        
        % implement concentrations for individual samples
        for a = 1:length(conc_individual_sample)
        
            % implement concentrations for all samples
            for i = 1:length(model.rxns)

                % if deltag data available
                if any(strcmp(data{1},model.rxns{i}))

                    % get concentrations of reactants and products
                    reactant_power = [];
                    reactant_lower = [];
                    reactant_upper = [];
                    product_power = [];
                    product_lower = [];
                    product_upper = [];

                    % iterate over reactants
                    for j = find(model.S(:,i) < 0)'
                        reactant_power = full(model.S(j,i));

                        % if metabolite found in concentration data
                        if any(strcmp(conc_individual_met{a},model.mets{j}))
                            reactant_lower(end+1) = conc_individual_lower{a}(strcmp(conc_individual_met{a},model.mets{j}));
                            reactant_upper(end+1) = conc_individual_upper{a}(strcmp(conc_individual_met{a},model.mets{j}));

                        % if metabolite not found in concentration data
                        else
                            reactant_lower(end+1) = default_concentration_lb;
                            reactant_upper(end+1) = default_concentration_ub;
                        end
                    end

                    % iterate over products
                    for j = find(model.S(:,i) > 0)'
                        product_power = full(model.S(j,i));

                        % if metabolite found in concentration data
                        if any(strcmp(conc_individual_met{a},model.mets{j}))
                            reactant_lower(end+1) = conc_individual_lower{a}(strcmp(conc_individual_met{a},model.mets{j}));
                            reactant_upper(end+1) = conc_individual_upper{a}(strcmp(conc_individual_met{a},model.mets{j}));

                        % if metabolite not found in concentration data
                        else
                            reactant_lower(end+1) = default_concentration_lb;
                            reactant_upper(end+1) = default_concentration_ub;
                        end
                    end

                    % calculate deltag bounds
                    Q_max = prod(product_upper.^product_power)/prod(reactant_lower.^reactant_power);
                    Q_min = prod(product_lower.^product_power)/prod(reactant_upper.^reactant_power);
                    deltag_lower(i,strcmp(samples_all,conc_individual_sample{a})) = data{4}(strcmp(data{1},model.rxns{i})) - 8.3144598/1000*310.15*log(Q_max)
                    deltag_upper(i,strcmp(samples_all,conc_individual_sample{a})) = data{4}(strcmp(data{1},model.rxns{i})) - 8.3144598/1000*310.15*log(Q_min)
                end
            end
        end
        
    % if not recalculating thermodynamics
    else
        
        % load deltaG file
        f = fopen('../data/deltag/deltag.csv','r');
        data = textscan(f,'%s %f %f %f %f %f','Delimiter',',','headerLines',1);
        fclose(f);
        
        % use provided minimum and maximum transformed deltag values
        deltag_lower = -Inf*ones(length(model.rxns),length(samples_all));
        deltag_upper = Inf*ones(length(model.rxns),length(samples_all));
        for i = 1:length(model.rxns)
            if any(strcmp(data{1},model.rxns{i}))
                deltag_lower(i,:) = data{5}(strcmp(data{1},model.rxns{i}));
                deltag_upper(i,:) = data{6}(strcmp(data{1},model.rxns{i}));
            end
        end
    end
end