
function [model,concentration_mets,concentration_values,constraints_fraction_id] = load_model(constraints_proliferation,constraints_proliferation_type,constraints_proliferation_value,constraints_value_type,constraints_value_value,constraints_value_reaction,constraints_value_count,constraints_fraction_type,constraints_fraction_value,constraints_fraction_reaction,constraints_fraction_count,constraints_maxmulti_weight,constraints_maxmulti_reaction,constraints_maxmulti_count,objective_function_direction,objective_function_weight,objective_function_reaction)

    % load model file
    load('../data/recon/recon3d_qflux.mat')
    
    % initialize concentrations
    concentration_mets = {};
    concentration_values = [];

    % compartment ph's
    concentration_mets = [concentration_mets,'h[c]','h[e]','h[g]','h[l]','h[m]','h[n]','h[r]','h[x]','h[i]'];
    concentration_values = [concentration_values,7.2,7.4,6.35,5.5,8,7.2,7.2,7,7];
    
    % implement model.c_fraction and model.c_obj to keep track of reactions that will have to maximize to add constraints
    model.c_fraction = zeros(size(model.c));
	model.c_maxmulti = zeros(size(model.c));
    model.c_obj = zeros(size(model.c));

    % implement custom value constraints
    for i = 1:constraints_value_count

        % if custom reaction
        if (contains(constraints_value_reaction{i},'-->')) || (contains(constraints_value_reaction{i},'<=>'))

            % initialize new reaction
            model.rxns{end+1} = sprintf('custom_value_%d',i);
            model.S(:,end+1) = zeros(length(model.mets),1);
            model.c(end+1) = 0;
            model.c_fraction(end+1) = 0;
			model.c_maxmulti(end+1) = 0;
            model.c_obj(end+1) = 0;
            model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
            model.rules{end+1} = '';
            model.grRules{end+1} = '';
            model.subSystems{end+1} = 'Custom Constraint';
            model.rxnNames{end+1} = sprintf('Custom Value Constraint %d',i);
            model.rxnKEGGID{end+1} = '';
            model.rxnKeggOrthology{end+1} = '';
            model.rxnConfidenceScores(end+1) = 0;
            model.rxnReferences{end+1} = '';
            model.rxnECNumbers{end+1} = '';
            model.rxnNotes{end+1} = '';
            model.rxnCOG{end+1} = '';
            model.rxnReconMap{end+1} = '';

            % irreversible reaction
            if (contains(constraints_value_reaction{i},'-->'))

                % parse reactants
                metstring = strsplit(strtrim(constraints_value_reaction{i}(1:(strfind(constraints_value_reaction{i},'-->')-1))),' + ','CollapseDelimiters',true);
                for j = 1:length(metstring)
                    if ~strcmp(metstring{j},'')
                        metstring_ind = strsplit(metstring{j},' ');

                        % metabolite
                        if ~any(strcmp(model.mets,metstring_ind{2}))
                            error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                        end

                        % stoichiometry
                        [num, status] = str2num(metstring_ind{1});
                        if status == 1
                            model.S(strcmp(model.mets,metstring_ind{2}),end) = -num;
                        else
                            error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                        end   
                    end
                end

                % parse products
                metstring = strsplit(strtrim(constraints_value_reaction{i}((strfind(constraints_value_reaction{i},'-->')+3):end)),' + ','CollapseDelimiters',true);
                for j = 1:length(metstring)
                    if ~strcmp(metstring{j},'')
                        metstring_ind = strsplit(metstring{j},' ');

                        % metabolite
                        if ~any(strcmp(model.mets,metstring_ind{2}))
                            error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                        end

                        % stoichiometry
                        [num, status] = str2num(metstring_ind{1});
                        if status == 1
                            model.S(strcmp(model.mets,metstring_ind{2}),end) = num;
                        else
                            error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                        end         
                    end
                end

                % check constraint value(s)
                for j = 1:length(constraints_value_value{i})
                    if constraints_value_value{i}(j) < 0
                        error('ERROR - Constraints - Custom - For irreversible reaction, constraint value must be >= 0')
                    end
                end

                % equal to
                if strcmp(constraints_value_type{i},'equal')
                    model.lb(end+1) = constraints_value_value{i}(1);
                    model.ub(end+1) = constraints_value_value{i}(1);

                % less than
                elseif strcmp(constraints_value_type{i},'less')
                    model.lb(end+1) = 0;
                    model.ub(end+1) = constraints_value_value{i}(1);

                % greater than
                elseif strcmp(constraints_value_type{i},'greater')
                    model.lb(end+1) = constraints_value_value{i}(1);
                    model.ub(end+1) = max(model.ub);

                % range
                elseif strcmp(constraints_value_type{i},'range')
                    if constraints_value_value{i}(1) <= constraints_value_value{i}(2)
                        model.lb(end+1) = constraints_value_value{i}(1);
                        model.ub(end+1) = constraints_value_value{i}(2);
                    else
                        error('ERROR - Constraints - Custom - For Range, min value must be <= max value');
                    end
                end
                
                % initialize formula
                model.rxnFormulas{end+1} = '';

                % substrates
                for j = find(model.S(:,end)<0)'
                    if isempty(model.rxnFormulas{end})
                        model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(j,end)),model.mets{j});
                    else
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(j,end)),model.mets{j})];
                    end
                end

                % reaction arrow
                if (model.lb(end) < 0) && (model.ub(end) > 0)
                    model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                elseif (model.lb(end) < 0)
                    model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                elseif (model.ub(end) > 0)
                    model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                else
                    model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                end

                % products
                numproducts = 0;
                for j = find(model.S(:,end)>0)'
                    if numproducts == 0
                        model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(j,end)),model.mets{j})];
                    else
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(j,end)),model.mets{j})];
                    end
                    numproducts = numproducts + 1;
                end

            % reversible reaction
            else

                % parse reactants
                metstring = strsplit(strtrim(constraints_value_reaction{i}(1:(strfind(constraints_value_reaction{i},'<=>')-1))),' + ','CollapseDelimiters',true);
                for j = 1:length(metstring)
                    if ~strcmp(metstring{j},'')
                        metstring_ind = strsplit(metstring{j},' ');

                        % metabolite
                        if ~any(strcmp(model.mets,metstring_ind{2}))
                            error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                        end

                        % stoichiometry
                        [num, status] = str2num(metstring_ind{1});
                        if status == 1
                            model.S(strcmp(model.mets,metstring_ind{2}),end) = -num;
                        else
                            error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                        end   
                    end
                end

                % parse products
                metstring = strsplit(strtrim(constraints_value_reaction{i}((strfind(constraints_value_reaction{i},'<=>')+3):end)),' + ','CollapseDelimiters',true);
                for j = 1:length(metstring)
                    if ~strcmp(metstring{j},'')
                        metstring_ind = strsplit(metstring{j},' ');

                        % metabolite
                        if ~any(strcmp(model.mets,metstring_ind{2}))
                            error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                        end

                        % stoichiometry
                        [num, status] = str2num(metstring_ind{1});
                        if status == 1
                            model.S(strcmp(model.mets,metstring_ind{2}),end) = num;
                        else
                            error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                        end       
                    end
                end

                % equal to
                if strcmp(constraints_value_type{i},'equal')
                    model.lb(end+1) = constraints_value_value{i}(1);
                    model.ub(end+1) = constraints_value_value{i}(1);

                % less than
                elseif strcmp(constraints_value_type{i},'less')
                    model.lb(end+1) = min(model.lb);
                    model.ub(end+1) = constraints_value_value{i}(1);

                % greater than
                elseif strcmp(constraints_value_type{i},'greater')
                    model.lb(end+1) = constraints_value_value{i}(1);
                    model.ub(end+1) = max(model.ub);

                % range
                elseif strcmp(constraints_value_type{i},'range')
                    if constraints_value_value{i}(1) <= constraints_value_value{i}(2)
                        model.lb(end+1) = constraints_value_value{i}(1);
                        model.ub(end+1) = constraints_value_value{i}(2);
                    else
                        error('ERROR - Constraints - Custom - For Range, min value must be <= max value');
                    end
                end
                
                % initialize formula
                model.rxnFormulas{end+1} = '';

                % substrates
                for j = find(model.S(:,end)<0)'
                    if isempty(model.rxnFormulas{end})
                        model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(j,end)),model.mets{j});
                    else
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(j,end)),model.mets{j})];
                    end
                end

                % reaction arrow
                if (model.lb(end) < 0) && (model.ub(end) > 0)
                    model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                elseif (model.lb(end) < 0)
                    model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                elseif (model.ub(end) > 0)
                    model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                else
                    model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                end

                % products
                numproducts = 0;
                for j = find(model.S(:,end)>0)'
                    if numproducts == 0
                        model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(j,end)),model.mets{j})];
                    else
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(j,end)),model.mets{j})];
                    end
                    numproducts = numproducts + 1;
                end
            end

        % if not custom reaction
        elseif any(strcmp(model.rxns,constraints_value_reaction{i}))

            % check if biomass_reaction
            if (strcmp(constraints_value_reaction{i},'biomass_reaction')) && (constraints_proliferation == 1)
                error('ERROR - Constraints - Custom - If already implementing proliferation constraint, cannot create additional constraints on biomass reaction');
            end

            % equal to
            if strcmp(constraints_value_type{i},'equal')
                model.lb(find(strcmp(model.rxns,constraints_value_reaction{i}))) = constraints_value_value{i}(1);
                model.ub(find(strcmp(model.rxns,constraints_value_reaction{i}))) = constraints_value_value{i}(1);

            % less than
            elseif strcmp(constraints_value_type{i},'less')
                model.ub(find(strcmp(model.rxns,constraints_value_reaction{i}))) = constraints_value_value{i}(1);

            % greater than
            elseif strcmp(constraints_value_type{i},'greater')
                model.lb(find(strcmp(model.rxns,constraints_value_reaction{i}))) = constraints_value_value{i}(1);

            % range
            elseif strcmp(constraints_value_type{i},'range')
                if constraints_value_value{i}(1) <= constraints_value_value{i}(2)
                    model.lb(find(strcmp(model.rxns,constraints_value_reaction{i}))) = constraints_value_value{i}(1);
                    model.ub(find(strcmp(model.rxns,constraints_value_reaction{i}))) = constraints_value_value{i}(2);
                else
                    error('ERROR - Constraints - Custom - For Range, min value must be <= max value');
                end
            end

        else
            error('ERROR - Constraints - Custom - Either invalid reaction formula or non-matching reaction name');
        end
    end
    
    % record rxn id's for custom fraction of maximum constraints
    constraints_fraction_id = [];
    
    % implement custom fraction of maximum constraints
    for i = 1:constraints_fraction_count
        
        % if custom reaction
        if (contains(constraints_fraction_reaction{i},'-->')) || (contains(constraints_fraction_reaction{i},'<=>'))

            % irreversible reaction
            if (contains(constraints_fraction_reaction{i},'-->'))

                % determine if all() reaction
                if strcmp(constraints_fraction_reaction{i}(1:4),'all(') && strcmp(constraints_fraction_reaction{i}(end),')')

                    % reaction formula without all()
                    rxn_no_all = constraints_fraction_reaction{i}(5:end-1);

                    % save reactant and product list
                    reactants = {};
                    reactant_stoich = [];
                    reactant_compartments = {};
                    products = {};
                    product_stoich = [];
                    product_compartments = {};

                    % parse reactants
                    metstring = strsplit(strtrim(rxn_no_all(1:(strfind(rxn_no_all,'-->')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                reactants{end+1} = metstring_ind{2}(1:end-2);
                                reactant_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    reactant_stoich(end+1) = num;
                                else
                                    error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Constraints - Custom - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(rxn_no_all((strfind(rxn_no_all,'-->')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                products{end+1} = metstring_ind{2}(1:end-2);
                                product_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    product_stoich(end+1) = num;
                                else
                                    error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Constraints - Custom - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % find compartments for all reactants and products
                    all_compartments = {'c','e','g','l','m','n','r','x','i'};
                    for j = 1:length(reactant_compartments)
                        all_compartments = intersect(all_compartments,reactant_compartments{j});
                    end
                    for j = 1:length(product_compartments)
                        all_compartments = intersect(all_compartments,product_compartments{j});
                    end

                    % if at least one compartment
                    if length(all_compartments) > 0

                        % create new reaction for all compartments
                        for j = 1:length(all_compartments)
                            model.rxns{end+1} = sprintf('custom_fraction_%d_%s',i,all_compartments{j});
							constraints_fraction_id(end+1) = length(model.rxns);
                            model.S(:,end+1) = zeros(length(model.mets),1);
                            model.lb(end+1) = 0;
                            model.ub(end+1) = max(model.ub);
                            model.c(end+1) = 0;
                            model.c_obj(end+1) = 0;
                            model.c_fraction(end+1) = i;
							model.c_maxmulti(end+1) = 0;
                            model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                            model.rules{end+1} = '';
                            model.grRules{end+1} = '';
                            model.subSystems{end+1} = 'Custom Constraint'
                            model.rxnNames{end+1} = sprintf('Custom Fraction Constraint %d - Compartment %s',i,all_compartments{j});
                            model.rxnKEGGID{end+1} = '';
                            model.rxnKeggOrthology{end+1} = '';
                            model.rxnConfidenceScores(end+1) = 0;
                            model.rxnReferences{end+1} = '';
                            model.rxnECNumbers{end+1} = '';
                            model.rxnNotes{end+1} = '';
                            model.rxnCOG{end+1} = '';
                            model.rxnReconMap{end+1} = '';

                            % reactants
                            for k = 1:length(reactants)
                                metstring = sprintf('%s[%s]',reactants{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = -reactant_stoich(k);
                            end

                            % products
                            for k = 1:length(products)
                                metstring = sprintf('%s[%s]',products{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = product_stoich(k);
                            end
                            
                            % initialize formula
                            model.rxnFormulas{end+1} = '';

                            % substrates
                            for k = find(model.S(:,end)<0)'
                                if isempty(model.rxnFormulas{end})
                                    model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                                end
                            end

                            % reaction arrow
                            if (model.lb(end) < 0) && (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                            elseif (model.lb(end) < 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                            elseif (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                            else
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                            end

                            % products
                            numproducts = 0;
                            for k = find(model.S(:,end)>0)'
                                if numproducts == 0
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                end
                                numproducts = numproducts + 1;
                            end
                        end    

                    % if no compartments
                    else
                       error('ERROR - Constraints - Custom - Could not find any compartments that contain all reaction metabolites') 
                    end

                % if not all() reaction
                else

                    % create new reaction
                    model.rxns{end+1} = sprintf('custom_fraction_%d',i);
					constraints_fraction_id(end+1) = length(model.rxns);
                    model.S(:,end+1) = zeros(length(model.mets),1);
                    model.lb(end+1) = 0;
                    model.ub(end+1) = max(model.ub);
                    model.c(end+1) = 0;
                    model.c_obj(end+1) = 0;
                    model.c_fraction(end+1) = i;
					model.c_maxmulti(end+1) = 0;
                    model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                    model.rules{end+1} = '';
                    model.grRules{end+1} = '';
                    model.subSystems{end+1} = 'Custom Constraint';
                    model.rxnNames{end+1} = sprintf('Custom Fraction Constraint %d',i);
                    model.rxnKEGGID{end+1} = '';
                    model.rxnKeggOrthology{end+1} = '';
                    model.rxnConfidenceScores(end+1) = 0;
                    model.rxnReferences{end+1} = '';
                    model.rxnECNumbers{end+1} = '';
                    model.rxnNotes{end+1} = '';
                    model.rxnCOG{end+1} = '';
                    model.rxnReconMap{end+1} = '';

                    % parse reactants
                    metstring = strsplit(strtrim(constraints_fraction_reaction{i}(1:(strfind(constraints_fraction_reaction{i},'-->')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = -num;
                            else
                                error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                            end  
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(constraints_fraction_reaction{i}((strfind(constraints_fraction_reaction{i},'-->')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = num;
                            else
                                error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                            end   
                        end
                    end
                    
                    % initialize formula
                    model.rxnFormulas{end+1} = '';

                    % substrates
                    for k = find(model.S(:,end)<0)'
                        if isempty(model.rxnFormulas{end})
                            model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                        end
                    end

                    % reaction arrow
                    if (model.lb(end) < 0) && (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                    elseif (model.lb(end) < 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                    elseif (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                    else
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                    end

                    % products
                    numproducts = 0;
                    for k = find(model.S(:,end)>0)'
                        if numproducts == 0
                            model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        end
                        numproducts = numproducts + 1;
                    end
                end

            % reversible reaction
            else
                % determine if all() reaction
                if strcmp(constraints_fraction_reaction{i}(1:4),'all(') && strcmp(constraints_fraction_reaction{i}(end),')')

                    % reaction formula without all()
                    rxn_no_all = constraints_fraction_reaction{i}(5:end-1);

                    % save reactant and product list
                    reactants = {};
                    reactant_stoich = [];
                    reactant_compartments = {};
                    products = {};
                    product_stoich = [];
                    product_compartments = {};

                    % parse reactants
                    metstring = strsplit(strtrim(rxn_no_all(1:(strfind(rxn_no_all,'<=>')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                reactants{end+1} = metstring_ind{2}(1:end-2);
                                reactant_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    reactant_stoich(end+1) = num;
                                else
                                    error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Constraints - Custom - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(rxn_no_all((strfind(rxn_no_all,'<=>')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                products{end+1} = metstring_ind{2}(1:end-2);
                                product_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    product_stoich(end+1) = num;
                                else
                                    error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Constraints - Custom - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % find compartments for all reactants and products
                    all_compartments = {'c','e','g','l','m','n','r','x','i'};
                    for j = 1:length(reactant_compartments)
                        all_compartments = intersect(all_compartments,reactant_compartments{j});
                    end
                    for j = 1:length(product_compartments)
                        all_compartments = intersect(all_compartments,product_compartments{j});
                    end

                    % if at least one compartment
                    if length(all_compartments) > 0

                        % create new reaction for all compartments
                        for j = 1:length(all_compartments)
                            model.rxns{end+1} = sprintf('custom_fraction_%d_%s',i,all_compartments{j});
							constraints_fraction_id(end+1) = length(model.rxns);
                            model.S(:,end+1) = zeros(length(model.mets),1);
                            model.lb(end+1) = min(model.lb);
                            model.ub(end+1) = max(model.ub);
                            model.c(end+1) = 0;
                            model.c_obj(end+1) = 0;
                            model.c_fraction(end+1) = i;
							model.c_maxmulti(end+1) = 0;
                            model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                            model.rules{end+1} = '';
                            model.grRules{end+1} = '';
                            model.subSystems{end+1} = 'Custom Constraint'
                            model.rxnNames{end+1} = sprintf('Custom Fraction Constraint %d - Compartment %s',i,all_compartments{j});
                            model.rxnKEGGID{end+1} = '';
                            model.rxnKeggOrthology{end+1} = '';
                            model.rxnConfidenceScores(end+1) = 0;
                            model.rxnReferences{end+1} = '';
                            model.rxnECNumbers{end+1} = '';
                            model.rxnNotes{end+1} = '';
                            model.rxnCOG{end+1} = '';
                            model.rxnReconMap{end+1} = '';

                            % reactants
                            for k = 1:length(reactants)
                                metstring = sprintf('%s[%s]',reactants{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = -reactant_stoich(k);
                            end

                            % products
                            for k = 1:length(products)
                                metstring = sprintf('%s[%s]',products{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = product_stoich(k);
                            end
                            
                            % initialize formula
                            model.rxnFormulas{end+1} = '';

                            % substrates
                            for k = find(model.S(:,end)<0)'
                                if isempty(model.rxnFormulas{end})
                                    model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                                end
                            end

                            % reaction arrow
                            if (model.lb(end) < 0) && (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                            elseif (model.lb(end) < 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                            elseif (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                            else
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                            end

                            % products
                            numproducts = 0;
                            for k = find(model.S(:,end)>0)'
                                if numproducts == 0
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                end
                                numproducts = numproducts + 1;
                            end
                        end    

                    % if no compartments
                    else
                       error('ERROR - Constraints - Custom - Could not find any compartments that contain all reaction metabolites') 
                    end

                % if not all() reaction
                else

                    % create new reaction
                    model.rxns{end+1} = sprintf('custom_fraction_%d',i);
					constraints_fraction_id(end+1) = length(model.rxns);
                    model.S(:,end+1) = zeros(length(model.mets),1);
                    model.lb(end+1) = min(model.lb);
                    model.ub(end+1) = max(model.ub);
                    model.c(end+1) = 0;
                    model.c_obj(end+1) = 0;
                    model.c_fraction(end+1) = i;
					model.c_maxmulti(end+1) = 0;
                    model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                    model.rules{end+1} = '';
                    model.grRules{end+1} = '';
                    model.subSystems{end+1} = 'Custom Constraint';
                    model.rxnNames{end+1} = sprintf('Custom Fraction Constraint %d',i);
                    model.rxnKEGGID{end+1} = '';
                    model.rxnKeggOrthology{end+1} = '';
                    model.rxnConfidenceScores(end+1) = 0;
                    model.rxnReferences{end+1} = '';
                    model.rxnECNumbers{end+1} = '';
                    model.rxnNotes{end+1} = '';
                    model.rxnCOG{end+1} = '';
                    model.rxnReconMap{end+1} = '';

                    % parse reactants
                    metstring = strsplit(strtrim(constraints_fraction_reaction{i}(1:(strfind(constraints_fraction_reaction{i},'<=>')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = -num;
                            else
                                error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                            end  
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(constraints_fraction_reaction{i}((strfind(constraints_fraction_reaction{i},'<=>')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = num;
                            else
                                error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                            end   
                        end
                    end
                    
                    % initialize formula
                    model.rxnFormulas{end+1} = '';

                    % substrates
                    for k = find(model.S(:,end)<0)'
                        if isempty(model.rxnFormulas{end})
                            model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                        end
                    end

                    % reaction arrow
                    if (model.lb(end) < 0) && (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                    elseif (model.lb(end) < 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                    elseif (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                    else
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                    end

                    % products
                    numproducts = 0;
                    for k = find(model.S(:,end)>0)'
                        if numproducts == 0
                            model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        end
                        numproducts = numproducts + 1;
                    end
                end
            end

        % if not custom reaction
        elseif any(strcmp(model.rxns,constraints_fraction_reaction{i}))
            constraints_fraction_id(end+1) = find(strcmp(model.rxns,constraints_fraction_reaction{i}));
            
            % check if biomass_reaction
            if (strcmp(constraints_fraction_reaction{i},'biomass_reaction')) && (constraints_proliferation == 1)
                error('ERROR - Constraints - Custom - If already implementing proliferation constraint, cannot create additional constraints on biomass reaction');
            end
            
            % change model.c_fraction
            model.c_fraction(strcmp(model.rxns,constraints_fraction_reaction{i})) = i;

        else
            error('ERROR - Constraints - Custom - Either invalid reaction formula or non-matching reaction name');
        end
 
    end

    % implement custom max multi constraints
    for i = 1:constraints_maxmulti_count
        
        % if custom reaction
        if (contains(constraints_maxmulti_reaction{i},'-->')) || (contains(constraints_maxmulti_reaction{i},'<=>'))

            % irreversible reaction
            if (contains(constraints_maxmulti_reaction{i},'-->'))

                % determine if all() reaction
                if strcmp(constraints_maxmulti_reaction{i}(1:4),'all(') && strcmp(constraints_maxmulti_reaction{i}(end),')')

                    % reaction formula without all()
                    rxn_no_all = constraints_maxmulti_reaction{i}(5:end-1);

                    % save reactant and product list
                    reactants = {};
                    reactant_stoich = [];
                    reactant_compartments = {};
                    products = {};
                    product_stoich = [];
                    product_compartments = {};

                    % parse reactants
                    metstring = strsplit(strtrim(rxn_no_all(1:(strfind(rxn_no_all,'-->')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                reactants{end+1} = metstring_ind{2}(1:end-2);
                                reactant_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    reactant_stoich(end+1) = num;
                                else
                                    error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Constraints - Custom - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(rxn_no_all((strfind(rxn_no_all,'-->')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                products{end+1} = metstring_ind{2}(1:end-2);
                                product_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    product_stoich(end+1) = num;
                                else
                                    error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Constraints - Custom - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % find compartments for all reactants and products
                    all_compartments = {'c','e','g','l','m','n','r','x','i'};
                    for j = 1:length(reactant_compartments)
                        all_compartments = intersect(all_compartments,reactant_compartments{j});
                    end
                    for j = 1:length(product_compartments)
                        all_compartments = intersect(all_compartments,product_compartments{j});
                    end

                    % if at least one compartment
                    if length(all_compartments) > 0

                        % create new reaction for all compartments
                        for j = 1:length(all_compartments)
                            model.rxns{end+1} = sprintf('custom_maxmulti_%d_%s',i,all_compartments{j});
                            model.S(:,end+1) = zeros(length(model.mets),1);
                            model.lb(end+1) = 0;
                            model.ub(end+1) = max(model.ub);
                            model.c(end+1) = 0;
                            model.c_obj(end+1) = 0;
							model.c_fraction(end+1) = 0;
                            model.c_maxmulti(end+1) = constraints_maxmulti_weight{i};
                            model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                            model.rules{end+1} = '';
                            model.grRules{end+1} = '';
                            model.subSystems{end+1} = 'Custom Constraint'
                            model.rxnNames{end+1} = sprintf('Custom Max Multi Constraint %d - Compartment %s',i,all_compartments{j});
                            model.rxnKEGGID{end+1} = '';
                            model.rxnKeggOrthology{end+1} = '';
                            model.rxnConfidenceScores(end+1) = 0;
                            model.rxnReferences{end+1} = '';
                            model.rxnECNumbers{end+1} = '';
                            model.rxnNotes{end+1} = '';
                            model.rxnCOG{end+1} = '';
                            model.rxnReconMap{end+1} = '';

                            % reactants
                            for k = 1:length(reactants)
                                metstring = sprintf('%s[%s]',reactants{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = -reactant_stoich(k);
                            end

                            % products
                            for k = 1:length(products)
                                metstring = sprintf('%s[%s]',products{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = product_stoich(k);
                            end
                            
                            % initialize formula
                            model.rxnFormulas{end+1} = '';

                            % substrates
                            for k = find(model.S(:,end)<0)'
                                if isempty(model.rxnFormulas{end})
                                    model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                                end
                            end

                            % reaction arrow
                            if (model.lb(end) < 0) && (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                            elseif (model.lb(end) < 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                            elseif (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                            else
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                            end

                            % products
                            numproducts = 0;
                            for k = find(model.S(:,end)>0)'
                                if numproducts == 0
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                end
                                numproducts = numproducts + 1;
                            end
                        end    

                    % if no compartments
                    else
                       error('ERROR - Constraints - Custom - Could not find any compartments that contain all reaction metabolites') 
                    end

                % if not all() reaction
                else

                    % create new reaction
                    model.rxns{end+1} = sprintf('custom_maxmulti_%d',i);
                    model.S(:,end+1) = zeros(length(model.mets),1);
                    model.lb(end+1) = 0;
                    model.ub(end+1) = max(model.ub);
                    model.c(end+1) = 0;
                    model.c_obj(end+1) = 0;
					model.c_fraction(end+1) = 0;
                    model.c_maxmulti(end+1) = constraints_maxmulti_weight{i};
                    model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                    model.rules{end+1} = '';
                    model.grRules{end+1} = '';
                    model.subSystems{end+1} = 'Custom Constraint';
                    model.rxnNames{end+1} = sprintf('Custom Max Multi Constraint %d',i);
                    model.rxnKEGGID{end+1} = '';
                    model.rxnKeggOrthology{end+1} = '';
                    model.rxnConfidenceScores(end+1) = 0;
                    model.rxnReferences{end+1} = '';
                    model.rxnECNumbers{end+1} = '';
                    model.rxnNotes{end+1} = '';
                    model.rxnCOG{end+1} = '';
                    model.rxnReconMap{end+1} = '';

                    % parse reactants
                    metstring = strsplit(strtrim(constraints_maxmulti_reaction{i}(1:(strfind(constraints_maxmulti_reaction{i},'-->')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = -num;
                            else
                                error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                            end  
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(constraints_maxmulti_reaction{i}((strfind(constraints_maxmulti_reaction{i},'-->')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = num;
                            else
                                error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                            end   
                        end
                    end
                    
                    % initialize formula
                    model.rxnFormulas{end+1} = '';

                    % substrates
                    for k = find(model.S(:,end)<0)'
                        if isempty(model.rxnFormulas{end})
                            model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                        end
                    end

                    % reaction arrow
                    if (model.lb(end) < 0) && (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                    elseif (model.lb(end) < 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                    elseif (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                    else
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                    end

                    % products
                    numproducts = 0;
                    for k = find(model.S(:,end)>0)'
                        if numproducts == 0
                            model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        end
                        numproducts = numproducts + 1;
                    end
                end

            % reversible reaction
            else
                % determine if all() reaction
                if strcmp(constraints_maxmulti_reaction{i}(1:4),'all(') && strcmp(constraints_maxmulti_reaction{i}(end),')')

                    % reaction formula without all()
                    rxn_no_all = constraints_maxmulti_reaction{i}(5:end-1);

                    % save reactant and product list
                    reactants = {};
                    reactant_stoich = [];
                    reactant_compartments = {};
                    products = {};
                    product_stoich = [];
                    product_compartments = {};

                    % parse reactants
                    metstring = strsplit(strtrim(rxn_no_all(1:(strfind(rxn_no_all,'<=>')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                reactants{end+1} = metstring_ind{2}(1:end-2);
                                reactant_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    reactant_stoich(end+1) = num;
                                else
                                    error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Constraints - Custom - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(rxn_no_all((strfind(rxn_no_all,'<=>')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                products{end+1} = metstring_ind{2}(1:end-2);
                                product_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    product_stoich(end+1) = num;
                                else
                                    error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Constraints - Custom - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % find compartments for all reactants and products
                    all_compartments = {'c','e','g','l','m','n','r','x','i'};
                    for j = 1:length(reactant_compartments)
                        all_compartments = intersect(all_compartments,reactant_compartments{j});
                    end
                    for j = 1:length(product_compartments)
                        all_compartments = intersect(all_compartments,product_compartments{j});
                    end

                    % if at least one compartment
                    if length(all_compartments) > 0

                        % create new reaction for all compartments
                        for j = 1:length(all_compartments)
                            model.rxns{end+1} = sprintf('custom_maxmulti_%d_%s',i,all_compartments{j});
                            model.S(:,end+1) = zeros(length(model.mets),1);
                            model.lb(end+1) = min(model.lb);
                            model.ub(end+1) = max(model.ub);
                            model.c(end+1) = 0;
                            model.c_obj(end+1) = 0;
							model.c_fraction(end+1) = 0;
                            model.c_maxmulti(end+1) = constraints_maxmulti_weight{i};
                            model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                            model.rules{end+1} = '';
                            model.grRules{end+1} = '';
                            model.subSystems{end+1} = 'Custom Constraint'
                            model.rxnNames{end+1} = sprintf('Custom Max Multi Constraint %d - Compartment %s',i,all_compartments{j});
                            model.rxnKEGGID{end+1} = '';
                            model.rxnKeggOrthology{end+1} = '';
                            model.rxnConfidenceScores(end+1) = 0;
                            model.rxnReferences{end+1} = '';
                            model.rxnECNumbers{end+1} = '';
                            model.rxnNotes{end+1} = '';
                            model.rxnCOG{end+1} = '';
                            model.rxnReconMap{end+1} = '';

                            % reactants
                            for k = 1:length(reactants)
                                metstring = sprintf('%s[%s]',reactants{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = -reactant_stoich(k);
                            end

                            % products
                            for k = 1:length(products)
                                metstring = sprintf('%s[%s]',products{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = product_stoich(k);
                            end
                            
                            % initialize formula
                            model.rxnFormulas{end+1} = '';

                            % substrates
                            for k = find(model.S(:,end)<0)'
                                if isempty(model.rxnFormulas{end})
                                    model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                                end
                            end

                            % reaction arrow
                            if (model.lb(end) < 0) && (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                            elseif (model.lb(end) < 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                            elseif (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                            else
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                            end

                            % products
                            numproducts = 0;
                            for k = find(model.S(:,end)>0)'
                                if numproducts == 0
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                end
                                numproducts = numproducts + 1;
                            end
                        end    

                    % if no compartments
                    else
                       error('ERROR - Constraints - Custom - Could not find any compartments that contain all reaction metabolites') 
                    end

                % if not all() reaction
                else

                    % create new reaction
                    model.rxns{end+1} = sprintf('custom_maxmulti_%d',i);
                    model.S(:,end+1) = zeros(length(model.mets),1);
                    model.lb(end+1) = min(model.lb);
                    model.ub(end+1) = max(model.ub);
                    model.c(end+1) = 0;
                    model.c_obj(end+1) = 0;
					model.c_fraction(end+1) = 0;
                    model.c_maxmulti(end+1) = constraints_maxmulti_weight{i};
                    model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                    model.rules{end+1} = '';
                    model.grRules{end+1} = '';
                    model.subSystems{end+1} = 'Custom Constraint';
                    model.rxnNames{end+1} = sprintf('Custom Max Multi Constraint %d',i);
                    model.rxnKEGGID{end+1} = '';
                    model.rxnKeggOrthology{end+1} = '';
                    model.rxnConfidenceScores(end+1) = 0;
                    model.rxnReferences{end+1} = '';
                    model.rxnECNumbers{end+1} = '';
                    model.rxnNotes{end+1} = '';
                    model.rxnCOG{end+1} = '';
                    model.rxnReconMap{end+1} = '';

                    % parse reactants
                    metstring = strsplit(strtrim(constraints_maxmulti_reaction{i}(1:(strfind(constraints_maxmulti_reaction{i},'<=>')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = -num;
                            else
                                error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                            end  
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(constraints_maxmulti_reaction{i}((strfind(constraints_maxmulti_reaction{i},'<=>')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Constraints - Custom - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = num;
                            else
                                error('ERROR - Constraints - Custom - Unrecognized stoichiometry %s',metstring_ind{1})
                            end   
                        end
                    end
                    
                    % initialize formula
                    model.rxnFormulas{end+1} = '';

                    % substrates
                    for k = find(model.S(:,end)<0)'
                        if isempty(model.rxnFormulas{end})
                            model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                        end
                    end

                    % reaction arrow
                    if (model.lb(end) < 0) && (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                    elseif (model.lb(end) < 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                    elseif (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                    else
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                    end

                    % products
                    numproducts = 0;
                    for k = find(model.S(:,end)>0)'
                        if numproducts == 0
                            model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        end
                        numproducts = numproducts + 1;
                    end
                end
            end

        % if not custom reaction
        elseif any(strcmp(model.rxns,constraints_maxmulti_reaction{i}))
            
            % change model.c_maxmulti
            model.c_maxmulti(strcmp(model.rxns,constraints_maxmulti_reaction{i})) = constraints_maxmulti_weight{i};

        else
            error('ERROR - Constraints - Custom - Either invalid reaction formula or non-matching reaction name');
        end
 
    end

    % implement objective functions
    for i = 1:length(objective_function_reaction)

        % if custom reaction
        if (contains(objective_function_reaction{i},'-->')) || (contains(objective_function_reaction{i},'<=>'))

            % irreversible reaction
            if (contains(objective_function_reaction{i},'-->'))

                % determine if all() reaction
                if strcmp(objective_function_reaction{i}(1:4),'all(') && strcmp(objective_function_reaction{i}(end),')')

                    % reaction formula without all()
                    rxn_no_all = objective_function_reaction{i}(5:end-1);

                    % save reactant and product list
                    reactants = {};
                    reactant_stoich = [];
                    reactant_compartments = {};
                    products = {};
                    product_stoich = [];
                    product_compartments = {};

                    % parse reactants
                    metstring = strsplit(strtrim(rxn_no_all(1:(strfind(rxn_no_all,'-->')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                reactants{end+1} = metstring_ind{2}(1:end-2);
                                reactant_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    reactant_stoich(end+1) = num;
                                else
                                    error('ERROR - Objective Function - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Objective Function - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(rxn_no_all((strfind(rxn_no_all,'-->')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                products{end+1} = metstring_ind{2}(1:end-2);
                                product_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    product_stoich(end+1) = num;
                                else
                                    error('ERROR - Objective Function - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Objective Function - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % find compartments for all reactants and products
                    all_compartments = {'c','e','g','l','m','n','r','x','i'};
                    for j = 1:length(reactant_compartments)
                        all_compartments = intersect(all_compartments,reactant_compartments{j});
                    end
                    for j = 1:length(product_compartments)
                        all_compartments = intersect(all_compartments,product_compartments{j});
                    end

                    % if at least one compartment
                    if length(all_compartments) > 0

                        % create new reaction for all compartments
                        for j = 1:length(all_compartments)
                            model.rxns{end+1} = sprintf('custom_objective_%d_%s',i,all_compartments{j});
                            model.S(:,end+1) = zeros(length(model.mets),1);
                            model.lb(end+1) = 0;
                            model.ub(end+1) = max(model.ub);
                            if strcmp(objective_function_direction{i},'MAX')
                                model.c(end+1) = objective_function_weight{i};
                                model.c_obj(end+1) = objective_function_weight{i};
                            elseif strcmp(objective_function_direction{i},'MIN')
                                model.c(end+1) = -objective_function_weight{i};
                                model.c_obj(end+1) = -objective_function_weight{i};
                            end
                            model.c_fraction(end+1) = 0;
							model.c_maxmulti(end+1) = 0;
                            model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                            model.rules{end+1} = '';
                            model.grRules{end+1} = '';
                            model.subSystems{end+1} = 'Custom Objective'
                            model.rxnNames{end+1} = sprintf('Custom Objective Function %d - Compartment %s',i,all_compartments{j});
                            model.rxnKEGGID{end+1} = '';
                            model.rxnKeggOrthology{end+1} = '';
                            model.rxnConfidenceScores(end+1) = 0;
                            model.rxnReferences{end+1} = '';
                            model.rxnECNumbers{end+1} = '';
                            model.rxnNotes{end+1} = '';
                            model.rxnCOG{end+1} = '';
                            model.rxnReconMap{end+1} = '';

                            % reactants
                            for k = 1:length(reactants)
                                metstring = sprintf('%s[%s]',reactants{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = -reactant_stoich(k);
                            end

                            % products
                            for k = 1:length(products)
                                metstring = sprintf('%s[%s]',products{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = product_stoich(k);
                            end
                            
                            % initialize formula
                            model.rxnFormulas{end+1} = '';

                            % substrates
                            for k = find(model.S(:,end)<0)'
                                if isempty(model.rxnFormulas{end})
                                    model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                                end
                            end

                            % reaction arrow
                            if (model.lb(end) < 0) && (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                            elseif (model.lb(end) < 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                            elseif (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                            else
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                            end

                            % products
                            numproducts = 0;
                            for k = find(model.S(:,end)>0)'
                                if numproducts == 0
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                end
                                numproducts = numproducts + 1;
                            end
                        end    

                    % if no compartments
                    else
                       error('ERROR - Objective Function - Could not find any compartments that contain all reaction metabolites') 
                    end

                % if not all() reaction
                else

                    % create new reaction
                    model.rxns{end+1} = sprintf('custom_objective_%d',i);
                    model.S(:,end+1) = zeros(length(model.mets),1);
                    model.lb(end+1) = 0;
                    model.ub(end+1) = max(model.ub);
                    if strcmp(objective_function_direction{i},'MAX')
                        model.c(end+1) = objective_function_weight{i};
                        model.c_obj(end+1) = objective_function_weight{i};
                    elseif strcmp(objective_function_direction{i},'MIN')
                        model.c(end+1) = -objective_function_weight{i};
                        model.c_obj(end+1) = -objective_function_weight{i};
                    end
                    model.c_fraction(end+1) = 0;
					model.c_maxmulti(end+1) = 0;
                    model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                    model.rules{end+1} = '';
                    model.grRules{end+1} = '';
                    model.subSystems{end+1} = 'Custom Objective';
                    model.rxnNames{end+1} = sprintf('Custom Objective Function %d',i);
                    model.rxnKEGGID{end+1} = '';
                    model.rxnKeggOrthology{end+1} = '';
                    model.rxnConfidenceScores(end+1) = 0;
                    model.rxnReferences{end+1} = '';
                    model.rxnECNumbers{end+1} = '';
                    model.rxnNotes{end+1} = '';
                    model.rxnCOG{end+1} = '';
                    model.rxnReconMap{end+1} = '';

                    % parse reactants
                    metstring = strsplit(strtrim(objective_function_reaction{i}(1:(strfind(objective_function_reaction{i},'-->')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Objective Function - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = -num;
                            else
                                error('ERROR - Objective Function - Unrecognized stoichiometry %s',metstring_ind{1})
                            end  
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(objective_function_reaction{i}((strfind(objective_function_reaction{i},'-->')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Objective Function - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = num;
                            else
                                error('ERROR - Objective Function - Unrecognized stoichiometry %s',metstring_ind{1})
                            end   
                        end
                    end
                    
                    % initialize formula
                    model.rxnFormulas{end+1} = '';

                    % substrates
                    for k = find(model.S(:,end)<0)'
                        if isempty(model.rxnFormulas{end})
                            model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                        end
                    end

                    % reaction arrow
                    if (model.lb(end) < 0) && (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                    elseif (model.lb(end) < 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                    elseif (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                    else
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                    end

                    % products
                    numproducts = 0;
                    for k = find(model.S(:,end)>0)'
                        if numproducts == 0
                            model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        end
                        numproducts = numproducts + 1;
                    end
                end

            % reversible reaction
            else
                % determine if all() reaction
                if strcmp(objective_function_reaction{i}(1:4),'all(') && strcmp(objective_function_reaction{i}(end),')')

                    % reaction formula without all()
                    rxn_no_all = objective_function_reaction{i}(5:end-1);

                    % save reactant and product list
                    reactants = {};
                    reactant_stoich = [];
                    reactant_compartments = {};
                    products = {};
                    product_stoich = [];
                    product_compartments = {};

                    % parse reactants
                    metstring = strsplit(strtrim(rxn_no_all(1:(strfind(rxn_no_all,'<=>')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                reactants{end+1} = metstring_ind{2}(1:end-2);
                                reactant_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    reactant_stoich(end+1) = num;
                                else
                                    error('ERROR - Objective Function - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Objective Function - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(rxn_no_all((strfind(rxn_no_all,'<=>')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % get matching metabolite and compartments
                            found_compartments = {};
                            for k = 1:length(model.mets)
                                if strcmp(model.mets{k}(1:end-3),metstring_ind{2}(1:end-2))
                                    found_compartments{end+1} = model.mets{k}(end-1);
                                end
                            end

                            % if found compartments
                            if length(found_compartments) > 0

                                % add metabolite
                                products{end+1} = metstring_ind{2}(1:end-2);
                                product_compartments{end+1} = found_compartments;

                                % add stoichiometry
                                [num, status] = str2num(metstring_ind{1});
                                if status == 1
                                    product_stoich(end+1) = num;
                                else
                                    error('ERROR - Objective Function - Unrecognized stoichiometry %s',metstring_ind{1})
                                end

                            % if no found compartments
                            else
                                error('ERROR - Objective Function - Could not find any compartments for metabolite %s',metstring_ind{1}(1:end-2))
                            end 
                        end
                    end

                    % find compartments for all reactants and products
                    all_compartments = {'c','e','g','l','m','n','r','x','i'};
                    for j = 1:length(reactant_compartments)
                        all_compartments = intersect(all_compartments,reactant_compartments{j});
                    end
                    for j = 1:length(product_compartments)
                        all_compartments = intersect(all_compartments,product_compartments{j});
                    end

                    % if at least one compartment
                    if length(all_compartments) > 0

                        % create new reaction for all compartments
                        for j = 1:length(all_compartments)
                            model.rxns{end+1} = sprintf('custom_objective_%d_%s',i,all_compartments{j});
                            model.S(:,end+1) = zeros(length(model.mets),1);
                            model.lb(end+1) = min(model.lb);
                            model.ub(end+1) = max(model.ub);
                            if strcmp(objective_function_direction{i},'MAX')
                                model.c(end+1) = objective_function_weight{i};
                                model.c_obj(end+1) = objective_function_weight{i};
                            elseif strcmp(objective_function_direction{i},'MIN')
                                model.c(end+1) = -objective_function_weight{i};
                                model.c_obj(end+1) = -objective_function_weight{i};
                            end
                            model.c_fraction(end+1) = 0;
							model.c_maxmulti(end+1) = 0;
                            model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                            model.rules{end+1} = '';
                            model.grRules{end+1} = '';
                            model.subSystems{end+1} = 'Custom Objective';
                            model.rxnNames{end+1} = sprintf('Custom Objective Function %d - Compartment %s',i,all_compartments{j});
                            model.rxnKEGGID{end+1} = '';
                            model.rxnKeggOrthology{end+1} = '';
                            model.rxnConfidenceScores(end+1) = 0;
                            model.rxnReferences{end+1} = '';
                            model.rxnECNumbers{end+1} = '';
                            model.rxnNotes{end+1} = '';
                            model.rxnCOG{end+1} = '';
                            model.rxnReconMap{end+1} = '';

                            % reactants
                            for k = 1:length(reactants)
                                metstring = sprintf('%s[%s]',reactants{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = -reactant_stoich(k);
                            end

                            % products
                            for k = 1:length(products)
                                metstring = sprintf('%s[%s]',products{k},all_compartments{j});
                                model.S(strcmp(model.mets,metstring),end) = product_stoich(k);
                            end
                            
                            % initialize formula
                            model.rxnFormulas{end+1} = '';

                            % substrates
                            for k = find(model.S(:,end)<0)'
                                if isempty(model.rxnFormulas{end})
                                    model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                                end
                            end

                            % reaction arrow
                            if (model.lb(end) < 0) && (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                            elseif (model.lb(end) < 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                            elseif (model.ub(end) > 0)
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                            else
                                model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                            end

                            % products
                            numproducts = 0;
                            for k = find(model.S(:,end)>0)'
                                if numproducts == 0
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                else
                                    model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                                end
                                numproducts = numproducts + 1;
                            end
                        end    

                    % if no compartments
                    else
                       error('ERROR - Objective Function - Could not find any compartments that contain all reaction metabolites') 
                    end

                % if not all() reaction
                else

                    % create new reaction
                    model.rxns{end+1} = sprintf('custom_objective_%d',i);
                    model.S(:,end+1) = zeros(length(model.mets),1);
                    model.lb(end+1) = min(model.lb);
                    model.ub(end+1) = max(model.ub);
                    if strcmp(objective_function_direction{i},'MAX')
                        model.c(end+1) = objective_function_weight{i};
                        model.c_obj(end+1) = objective_function_weight{i};
                    elseif strcmp(objective_function_direction{i},'MIN')
                        model.c(end+1) = -objective_function_weight{i};
                        model.c_obj(end+1) = -objective_function_weight{i};
                    end
                    model.c_fraction(end+1) = 0;
					model.c_maxmulti(end+1) = 0;
                    model.rxnGeneMat(end+1,:) = zeros(1,length(model.genes));
                    model.rules{end+1} = '';
                    model.grRules{end+1} = '';
                    model.subSystems{end+1} = 'Custom Objective';
                    model.rxnNames{end+1} = 'Custom Objective Function';
                    model.rxnKEGGID{end+1} = '';
                    model.rxnKeggOrthology{end+1} = '';
                    model.rxnConfidenceScores(end+1) = 0;
                    model.rxnReferences{end+1} = '';
                    model.rxnECNumbers{end+1} = '';
                    model.rxnNotes{end+1} = '';
                    model.rxnCOG{end+1} = '';
                    model.rxnReconMap{end+1} = '';

                    % parse reactants
                    metstring = strsplit(strtrim(objective_function_reaction{i}(1:(strfind(objective_function_reaction{i},'<=>')-1))),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Objective Function - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = -num;
                            else
                                error('ERROR - Objective Function - Unrecognized stoichiometry %s',metstring_ind{1})
                            end  
                        end
                    end

                    % parse products
                    metstring = strsplit(strtrim(objective_function_reaction{i}((strfind(objective_function_reaction{i},'<=>')+3):end)),' + ','CollapseDelimiters',true);
                    for j = 1:length(metstring)
                        if ~strcmp(metstring{j},'')
                            metstring_ind = strsplit(metstring{j},' ');

                            % metabolite
                            if ~any(strcmp(model.mets,metstring_ind{2}))
                                error('ERROR - Objective Function - Unrecognized metabolite %s',metstring_ind{2})
                            end

                            % stoichiometry
                            [num, status] = str2num(metstring_ind{1});
                            if status == 1
                                model.S(strcmp(model.mets,metstring_ind{2}),end) = num;
                            else
                                error('ERROR - Objective Function - Unrecognized stoichiometry %s',metstring_ind{1})
                            end   
                        end
                    end
                    
                    % initialize formula
                    model.rxnFormulas{end+1} = '';

                    % substrates
                    for k = find(model.S(:,end)<0)'
                        if isempty(model.rxnFormulas{end})
                            model.rxnFormulas{end} = sprintf('%d %s',-full(model.S(k,end)),model.mets{k});
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',-full(model.S(k,end)),model.mets{k})];
                        end
                    end

                    % reaction arrow
                    if (model.lb(end) < 0) && (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <=> '];
                    elseif (model.lb(end) < 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' <-- '];
                    elseif (model.ub(end) > 0)
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' --> '];
                    else
                        model.rxnFormulas{end} = [model.rxnFormulas{end},' X--X '];
                    end

                    % products
                    numproducts = 0;
                    for k = find(model.S(:,end)>0)'
                        if numproducts == 0
                            model.rxnFormulas{end} = [model.rxnFormulas{end},sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        else
                            model.rxnFormulas{end} = [model.rxnFormulas{end},' + ',sprintf('%d %s',full(model.S(k,end)),model.mets{k})];
                        end
                        numproducts = numproducts + 1;
                    end
                end
            end

        % if not custom reaction
        elseif any(strcmp(model.rxns,objective_function_reaction{i}))

            if strcmp(objective_function_direction{i},'MAX')
                model.c(strcmp(model.rxns,objective_function_reaction{i})) = objective_function_weight{i};
                model.c_obj(strcmp(model.rxns,objective_function_reaction{i})) = objective_function_weight{i};
            elseif strcmp(objective_function_direction{i},'MIN')
                model.c(strcmp(model.rxns,objective_function_reaction{i})) = -objective_function_weight{i};
                model.c_obj(strcmp(model.rxns,objective_function_reaction{i})) = -objective_function_weight{i};
            end

        else
            error('ERROR - Objective Function - Either invalid reaction formula or non-matching reaction name');
        end

    end
end