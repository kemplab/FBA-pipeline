
function [] = run_fba(output_folder,input_file,model,type_of_analysis,samples_all,samples_all_source,samples_all_source_protein,vmax_forward,vmax_reverse,objective_function_direction,objective_function_weight,objective_function_reaction,exchange_rxn,deltag_lower,deltag_upper,constraints_fraction_count,constraints_fraction_type,constraints_fraction_reaction,constraints_fraction_value,constraints_fraction_id,constraints_maxmulti_fraction,constraints_maxmulti_count,pfba_fraction,fva_all_fraction,fva_select_fraction,fva_select_groups,fva_select_reactions,sampler_fraction,sampler_samples,sampler_steps,media_sensitivity_critical_fraction,media_sensitivity_critical_metabolite1,media_sensitivity_spectrum_fraction,media_sensitivity_spectrum_metabolite1,media_sensitivity_spectrum_metabolite2,media_sensitivity_spectrum_points,parallel_cores,knockdown_genes,knockdown_fractions,catalase_version,module_5fu,module_cis,module_cpa,module_dox);

    % initialize results
    temp = strsplit(input_file,'.xlsx');
    temp = strsplit(temp{1},'/');
    input_filename = temp{end};
    mkdir(sprintf('results/%s/%s',output_folder,input_filename))
    if type_of_analysis == 2
        mkdir(sprintf('results/%s/%s/pfba',output_folder,input_filename))
    elseif (type_of_analysis == 3) || (type_of_analysis == 4)
        mkdir(sprintf('results/%s/%s/fva',output_folder,input_filename))
    elseif type_of_analysis == 5
        mkdir(sprintf('results/%s/%s/sampling',output_folder,input_filename))
    elseif type_of_analysis == 7
        mkdir(sprintf('results/%s/%s/media_sensitivity',output_folder,input_filename))
    end
    
    % create objective value file
    f_obj = fopen(sprintf('results/%s/%s/objval.tsv',output_folder,input_filename),'w');
    fprintf(f_obj,'SAMPLE\tSOURCE\tOBJVAL\n');
        
    % create critical media sensitivity file
    if type_of_analysis == 6
        f_media = fopen(sprintf('results/%s/%s/media_sensitivity_critical.tsv',output_folder,input_filename),'w');
        fprintf(f_media,'SAMPLE\tSOURCE\tCRITICAL UPTAKE [mmol/gDW/hr]\n');
    end

    % save initial model object
    modelOriginal = model;
    
    % iterate over each sample
    for i = 1:length(samples_all)
        
        % load original model
        model = modelOriginal;

        % implement vmax values
        for j = 1:length(model.rxns)

            % forward value
            if (model.ub(j) == max(model.ub)) && (~(isnan(vmax_forward(j,i))))
                model.ub(j) = vmax_forward(j,i);
            end

            % reverse value
            if (model.lb(j) == min(model.lb)) && (~(isnan(vmax_reverse(j,i))))
                model.lb(j) = -vmax_reverse(j,i);
            end   
        end

        % implement deltag constraints
        for j = 1:length(model.rxns)

            % if deltag range only positive, restrict forward direction
            if deltag_lower(j,i) > 0
                model.ub(j) = 0;
            end

            % if deltag range only negative, restrict reverse direction (if reaction is reversible)
            if deltag_upper(j,i) < 0
                model.lb(j) = 0;
            end
        end
        
        % implement catalase
        catalase_reactions = {'r0010','CATm','CATp','CATr'};
        if catalase_version == 0
            for j = 1:length(catalase_reactions)
                model.lb(strcmp(model.rxns,strcat(catalase_reactions{j},'_nadph'))) = 0;
                model.ub(strcmp(model.rxns,strcat(catalase_reactions{j},'_nadph'))) = 0;
            end
        else
            for j = 1:length(catalase_reactions)
                model.lb(strcmp(model.rxns,catalase_reactions{j})) = 0;
                model.ub(strcmp(model.rxns,catalase_reactions{j})) = 0;
            end
        end
        
        % implement 5-fluorouracil module
        if module_5fu == 0
            for j = 1:length(model.rxns)
                if strcmp(model.subSystems{j},'MODULE 5-FLUOROURACIL');
                    model.lb(j) = 0;
                    model.ub(j) = 0;
                end
            end
        end
        
        % implement cisplatin module
        if module_cis == 0
            for j = 1:length(model.rxns)
                if strcmp(model.subSystems{j},'MODULE CISPLATIN');
                    model.lb(j) = 0;
                    model.ub(j) = 0;
                end
            end
        end
        
        % implement cyclophosphamide module
        if module_cpa == 0
            for j = 1:length(model.rxns)
                if strcmp(model.subSystems{j},'MODULE CYCLOPHOSPHAMIDE');
                    model.lb(j) = 0;
                    model.ub(j) = 0;
                end
            end
        end
        
        % implement doxorubicin module
        if module_dox == 0
            for j = 1:length(model.rxns)
                if strcmp(model.subSystems{j},'MODULE DOXORUBICIN');
                    model.lb(j) = 0;
                    model.ub(j) = 0;
                end
            end
        end
        
        % setup LP problem
        model.A = model.S;
        model.rhs = model.b;
        model.obj = model.c;
        model.sense = repmat('=',1,length(model.mets));
        model.vtype = repmat('C',1,length(model.rxns));
        model.varnames = model.rxns;
        model.modelsense = lower(objective_function_direction);

        % reset objective value
        model.c = zeros(size(model.c));
        model.obj = zeros(size(model.obj));

        % find fraction of maximum constraints
        max_constraint_values = [];
        for j = 1:constraints_fraction_count

            % implement objective value
            model.c(model.c_fraction == j) = 1;
            model.obj(model.c_fraction == j) = 1;
            model.modelsense = 'max';

            % get objective value
            params.outputflag = 0;
            result = gurobi(model,params);
            max_constraint_values(end+1) = result.objval;

            % reset objective value
            model.c = zeros(size(model.c));
            model.obj = zeros(size(model.obj));
        end

        % implement fraction of maximum constraints
        for j = 1:constraints_fraction_count

            % equal to
            if strcmp(constraints_fraction_type{j},'equal')
                model.lb(constraints_fraction_id(j)) = constraints_fraction_value{j}(1) * max_constraint_values(j);
                model.ub(constraints_fraction_id(j)) = constraints_fraction_value{j}(1) * max_constraint_values(j);

            % less than
            elseif strcmp(constraints_fraction_type{j},'less')
                model.ub(constraints_fraction_id(j)) = constraints_fraction_value{j}(1) * max_constraint_values(j);

            % greater than
            elseif strcmp(constraints_fraction_type{j},'greater')
                model.lb(constraints_fraction_id(j)) = constraints_fraction_value{j}(1) * max_constraint_values(j);

            % range
            elseif strcmp(constraints_fraction_type{j},'range')
                if constraints_fraction_value{j}(1) <= constraints_fraction_value{j}(2)
                    model.lb(constraints_fraction_id(j)) = constraints_fraction_value{j}(1) * max_constraint_values(j);
                    model.ub(constraints_fraction_id(j)) = constraints_fraction_value{j}(2) * max_constraint_values(j);
                else
                    error('ERROR - Constraints - Custom - For Range, min value must be <= max value');
                end
            end      
        end  

		% find max multi constraints
		if constraints_maxmulti_count > 0
		
			% implement objective value
			model.c = model.c_maxmulti;
			model.obj = model.c_maxmulti;
			model.modelsense = 'max';

            % get objective value
            params.outputflag = 0;
            result = gurobi(model,params);
            maxmulti_value = result.objval;

            % reset objective value
            model.c = zeros(size(model.c));
            model.obj = zeros(size(model.obj));
        end

        % implement objective function
        model.c = model.c_obj;
        model.obj = model.c_obj;
        model.modelsense = 'max';
		
		% implement max multi constraint
		if constraints_maxmulti_count > 0
			model.S(end+1,:) = model.c_maxmulti;
			model.b(end+1) = constraints_maxmulti_fraction*maxmulti_value;
			model.csense(end+1) = 'E';
			model.mets{end+1} = 'maxmulti';
			model.metCharges(end+1) = 0;
			model.metFormulas{end+1} = '';
			model.metSmiles{end+1} = '';
			model.metNames{end+1} = 'Max Multi Constraint';
			model.metHMDBID{end+1} = '';
			model.metInChIString{end+1} = '';
			model.metKEGGID{end+1} = '';
			model.metPubChemID{end+1} = '';
			model.metCHEBIID{end+1} = '';
			model.metPdMap{end+1} = '';
			model.metReconMap{end+1} = '';
		end
		
		% create parameters
		model.A = model.S;
		model.rhs = model.b;
		model.sense = repmat('=',1,length(model.csense));
		for j = 1:length(model.csense)
			if model.csense(i) == 'G'
				model.sense(i) = '>';
			elseif model.csense(i) == 'L'
				model.sense(i) = '<';
			end
		end
		
        % objective function value
        params.outputflag = 0;
        result = gurobi(model,params);
        fprintf(f_obj,'%s\t%s\t%0.9f\n',samples_all{i},samples_all_source_protein{i},result.objval);

        % parsimonious fba
        if type_of_analysis == 2

            % ensure only one maximization objective
            if length(objective_function_reaction) ~= 1
                error('ERROR - General - Type of Analysis - Can only perform pFBA with maximizing one objective function reaction');
            elseif ~strcmp(objective_function_direction{1},'MAX')
                error('ERROR - General - Type of Analysis - Can only perform pFBA with maximizing one objective function reaction');
            end

            % create irreversible model
            model_irrev = model;
            for j = 1:length(model_irrev.rxns)

                % if both directions, split into two separate reactions
                if (model_irrev.lb(j) < 0) && (model_irrev.ub(j) > 0)

                    % create new reaction
                    model_irrev.rxns{end+1} = sprintf('%s_REV',model_irrev.rxns{j});
                    model_irrev.S(:,end+1) = -model_irrev.S(:,j);
                    model_irrev.lb(end+1) = 0;
                    model_irrev.ub(end+1) = -model_irrev.lb(j);
                    model_irrev.c(end+1) = 0;
                    model_irrev.rxnGeneMat(end+1,:) = model_irrev.rxnGeneMat(j,:);
                    model_irrev.rules{end+1} = model_irrev.rules{j};
                    model_irrev.grRules{end+1} = model_irrev.grRules{j};
                    model_irrev.subSystems{end+1} = model_irrev.subSystems{j};
                    model_irrev.rxnNames{end+1} = model_irrev.rxnNames{j};
                    model_irrev.rxnKEGGID{end+1} = model_irrev.rxnKEGGID{j};
                    model_irrev.rxnKeggOrthology{end+1} = model_irrev.rxnKeggOrthology{j};
                    model_irrev.rxnConfidenceScores(end+1) = model_irrev.rxnConfidenceScores(j);
                    model_irrev.rxnReferences{end+1} = model_irrev.rxnReferences{j};
                    model_irrev.rxnECNumbers{end+1} = model_irrev.rxnECNumbers{j};
                    model_irrev.rxnNotes{end+1} = model_irrev.rxnNotes{j};
                    model_irrev.rxnCOG{end+1} = model_irrev.rxnCOG{j};
                    model_irrev.rxnReconMap{end+1} = model_irrev.rxnReconMap{j};

                    % make original reaction irreversible
                    model_irrev.lb(j) = 0;                

                % if just reverse direction, reverse to forward direction
                elseif (model_irrev.lb(j) < 0)

                    % reverse reaction
                    model_irrev.rxns{j} = sprintf('%s_REV',model_irrev.rxns{j});
                    model_irrev.S(:,j) = -model_irrev.S(:,j);
                    model_irrev.ub(j) = -model_irrev.lb(j);
                    model_irrev.lb(j) = 0;

                end
            end

            % implement objective value constraint
            model_irrev.S(end+1,:) = zeros(1,length(model_irrev.rxns));
            model_irrev.S(end,find(model_irrev.c)) = 1;
            model_irrev.b(end+1) = pfba_fraction * result.objval;

            % implement minimization of total flux
            original_model_irrev_c = model_irrev.c;
            model_irrev.c = zeros(size(model_irrev.c));
            model_irrev.c(~original_model_irrev_c) = 1;

            % setup LP problem
            model_irrev.A = model_irrev.S;
            model_irrev.rhs = model_irrev.b;
            model_irrev.obj = model_irrev.c;
            model_irrev.sense = repmat('=',1,length(model_irrev.mets)+1);
            model_irrev.vtype = repmat('C',1,length(model_irrev.rxns));
            model_irrev.varnames = model_irrev.rxns;
            model_irrev.modelsense = 'min';

            % pFBA fluxes
            params.outputflag = 0;
            result = gurobi(model_irrev,params);

            % collect flux values
            flux_value = [];
            flux_lb = [];
            flux_ub = [];
            for j = 1:length(model.rxns)
                if (any(strcmp(model_irrev.rxns,model.rxns{j}))) && (any(strcmp(model_irrev.rxns,sprintf('%s_REV',model.rxns{j}))))
                    if (result.x(strcmp(model_irrev.rxns,model.rxns{j})) > 0) && (result.x(strcmp(model_irrev.rxns,sprintf('%s_REV',model.rxns{j}))) > 0)
                        error(sprintf('pFBA - both directions of %s have a nonzero flux',model.rxns{j}));
                    elseif result.x(strcmp(model_irrev.rxns,sprintf('%s_REV',model.rxns{j}))) > 0
                        flux_value(end+1) = -result.x(strcmp(model_irrev.rxns,sprintf('%s_REV',model.rxns{j})));
                    else
                        flux_value(end+1) = result.x(strcmp(model_irrev.rxns,model.rxns{j}));
                    end
                    flux_lb(end+1) = -model_irrev.ub(strcmp(model_irrev.rxns,sprintf('%s_REV',model.rxns{j})));
                    flux_ub(end+1) = model_irrev.ub(strcmp(model_irrev.rxns,model.rxns{j}));
                elseif any(strcmp(model_irrev.rxns,sprintf('%s_REV',model.rxns{j})))
                    flux_value(end+1) = -result.x(strcmp(model_irrev.rxns,sprintf('%s_REV',model.rxns{j})));
                    flux_lb(end+1) = -model_irrev.ub(strcmp(model_irrev.rxns,sprintf('%s_REV',model.rxns{j})));
                    flux_ub(end+1) = 0;
                else
                    flux_value(end+1) = result.x(strcmp(model_irrev.rxns,model.rxns{j}));
                    flux_lb(end+1) = 0;
                    flux_ub(end+1) = model_irrev.ub(strcmp(model_irrev.rxns,model.rxns{j}));
                end
            end

            % output results
            f_pfba = fopen(sprintf('results/%s/%s/pfba/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
            fprintf(f_pfba,'REACTION\tFLUX [mmol/gDW/hr]\tLB [mmol/gDW/hr]\tUB [mmol/gDW/hr]\tNAME\tFORMULA\tSUBSYSTEM\tGR RULE\tEC NUMBER\tKEGG\tNOTES\tREFERENCES\n');
            for j = 1:length(model.rxns)
                fprintf(f_pfba,'%s\t%0.9f\t%0.9f\t%0.9f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',model.rxns{j},flux_value(j),flux_lb(j),flux_ub(j),model.rxnNames{j},model.rxnFormulas{j},model.subSystems{j},model.grRules{j},model.rxnECNumbers{j},model.rxnKEGGID{j},model.rxnNotes{j},model.rxnReferences{j});
            end
            fclose(f_pfba);

        % flux variability analysis - all reactions
        elseif type_of_analysis == 3

            % ensure only one maximization objective
            if length(objective_function_reaction) ~= 1
                error('ERROR - General - Type of Analysis - Can only perform FVA with maximizing one objective function reaction');
            elseif ~strcmp(objective_function_direction{1},'MAX')
                error('ERROR - General - Type of Analysis - Can only perform FVA with maximizing one objective function reaction');
            end
            
            % implement objective value constraint
            model.S(end+1,:) = zeros(1,length(model.rxns));
            model.S(end,find(model.c)) = 1;
            model.b(end+1) = fva_all_fraction * result.objval;

            % clear model.c
            model.c = zeros(size(model.c));

            % calculate FVA fluxes
            flux_min = [];
            flux_max = [];
            for j = 1:length(model.rxns)

                % set model.c
                model.c(j) = 1;

                % min value
                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = [repmat('=',1,length(model.mets)),'>'];
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'min';

                params.outputflag = 0;
                result = gurobi(model,params);
                flux_min(end+1) = result.objval;

                % max value
                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = [repmat('=',1,length(model.mets)),'>'];
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                flux_max(end+1) = result.objval;

                % clear model.c
                model.c = zeros(size(model.c));
            end

            % output results
            f_fva = fopen(sprintf('results/%s/%s/fva/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
            fprintf(f_fva,'REACTION\tMIN [mmol/gDW/hr]\tMAX [mmol/gDW/hr]\tLB [mmol/gDW/hr]\tUB [mmol/gDW/hr]\tNAME\tFORMULA\tSUBSYSTEM\tGR RULE\tEC NUMBER\tKEGG\tNOTES\tREFERENCES\n');
            for j = 1:length(model.rxns)
                fprintf(f_fva,'%s\t%0.9f\t%0.9f\t%0.9f\t%0.9f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',model.rxns{j},flux_min(j),flux_max(j),model.lb(j),model.ub(j),model.rxnNames{j},model.rxnFormulas{j},model.subSystems{j},model.grRules{j},model.rxnECNumbers{j},model.rxnKEGGID{j},model.rxnNotes{j},model.rxnReferences{j});
            end
            fclose(f_fva);

        % flux variability analysis - select reactions
        elseif type_of_analysis == 4

            % ensure only one maximization objective
            if length(objective_function_reaction) ~= 1
                error('ERROR - General - Type of Analysis - Can only perform FVA with maximizing one objective function reaction');
            elseif ~strcmp(objective_function_direction{1},'MAX')
                error('ERROR - General - Type of Analysis - Can only perform FVA with maximizing one objective function reaction');
            end

            % reaction list
            rxnlist = [];
            for j = 1:length(fva_select_reactions)
                if any(strcmp(model.rxns,fva_select_reactions{j}))
                    rxnlist(end+1) = find(strcmp(model.rxns,fva_select_reactions{j}));
                else
                    error('ERROR - General - Options- FVA - Select reactions - Reaction "%s" not found in model',fva_select_reactions{j});
                end
            end

            % add objective reactions
            for j = 1:length(model.rxns)
                if model.c(j) && ~any(j == rxnlist)
                    rxnlist(end+1) = j;
                end
            end

            % reaction groups
            for j = 1:length(fva_select_groups)
                switch fva_select_groups(j)

                    % all nadp+ --> nadph reactions
                    case 1

                        % find nadp
                        nadp = [];
                        for k = 1:length(model.mets)
                            if length(model.mets{k}) == 7
                                if strcmp(model.mets{k}(1:5),'nadp[')
                                    nadp(end+1) = k;
                                end
                            end
                        end

                        % find nadph
                        nadph = [];
                        for k = 1:length(model.mets)
                            if length(model.mets{k}) == 8
                                if strcmp(model.mets{k}(1:6),'nadph[')
                                    nadph(end+1) = k;
                                end
                            end
                        end

                        % find reactions
                        for k = 1:length(model.rxns)

                            % check forward direction
                            if model.ub(k) > 0

                                % find reactant
                                reactant_found = false;
                                for m = find(full(model.S(:,k)) < 0)'
                                    if any(m == nadp)
                                        reactant_found = true;
                                        break;
                                    end
                                end

                                % find product
                                product_found = false;
                                for m = find(full(model.S(:,k)) > 0)'
                                    if any(m == nadph)
                                        product_found = true;
                                        break;
                                    end
                                end

                                % add
                                if reactant_found && product_found && ~any(k == rxnlist)
                                    rxnlist(end+1) = k;
                                end
                            end

                            % check reverse direction
                            if model.ub(k) < 0

                                % find reactant
                                reactant_found = false;
                                for m = find(full(model.S(:,k)) > 0)'
                                    if any(m == nadp)
                                        reactant_found = true;
                                        break;
                                    end
                                end

                                % find product
                                product_found = false;
                                for m = find(full(model.S(:,k)) < 0)'
                                    if any(m == nadph)
                                        product_found = true;
                                        break;
                                    end
                                end

                                % add
                                if reactant_found && product_found && ~any(k == rxnlist)
                                    rxnlist(end+1) = k;
                                end
                            end
                        end


                    % all nadph --> nadp+ reactions
                    case 2

                        % find nadp
                        nadp = [];
                        for k = 1:length(model.mets)
                            if length(model.mets{k}) == 7
                                if strcmp(model.mets{k}(1:5),'nadp[')
                                    nadp(end+1) = k;
                                end
                            end
                        end

                        % find nadph
                        nadph = [];
                        for k = 1:length(model.mets)
                            if length(model.mets{k}) == 8
                                if strcmp(model.mets{k}(1:6),'nadph[')
                                    nadph(end+1) = k;
                                end
                            end
                        end

                        % find reactions
                        for k = 1:length(model.rxns)

                            % check forward direction
                            if model.ub(k) > 0

                                % find reactant
                                reactant_found = false;
                                for m = find(full(model.S(:,k)) < 0)'
                                    if any(m == nadph)
                                        reactant_found = true;
                                        break;
                                    end
                                end

                                % find product
                                product_found = false;
                                for m = find(full(model.S(:,k)) > 0)'
                                    if any(m == nadp)
                                        product_found = true;
                                        break;
                                    end
                                end

                                % add
                                if reactant_found && product_found && ~any(k == rxnlist)
                                    rxnlist(end+1) = k;
                                end
                            end

                            % check reverse direction
                            if model.ub(k) < 0

                                % find reactant
                                reactant_found = false;
                                for m = find(full(model.S(:,k)) > 0)'
                                    if any(m == nadph)
                                        reactant_found = true;
                                        break;
                                    end
                                end

                                % find product
                                product_found = false;
                                for m = find(full(model.S(:,k)) < 0)'
                                    if any(m == nadp)
                                        product_found = true;
                                        break;
                                    end
                                end

                                % add
                                if reactant_found && product_found && ~any(k == rxnlist)
                                    rxnlist(end+1) = k;
                                end
                            end
                        end

                    % all nad+ --> nadh reactions
                    case 3

                        % find nad
                        nad = [];
                        for k = 1:length(model.mets)
                            if length(model.mets{k}) == 6
                                if strcmp(model.mets{k}(1:4),'nad[')
                                    nad(end+1) = k;
                                end
                            end
                        end

                        % find nadh
                        nadh = [];
                        for k = 1:length(model.mets)
                            if length(model.mets{k}) == 7
                                if strcmp(model.mets{k}(1:5),'nadh[')
                                    nadh(end+1) = k;
                                end
                            end
                        end

                        % find reactions
                        for k = 1:length(model.rxns)

                            % check forward direction
                            if model.ub(k) > 0

                                % find reactant
                                reactant_found = false;
                                for m = find(full(model.S(:,k)) < 0)'
                                    if any(m == nad)
                                        reactant_found = true;
                                        break;
                                    end
                                end

                                % find product
                                product_found = false;
                                for m = find(full(model.S(:,k)) > 0)'
                                    if any(m == nadh)
                                        product_found = true;
                                        break;
                                    end
                                end

                                % add
                                if reactant_found && product_found && ~any(k == rxnlist)
                                    rxnlist(end+1) = k;
                                end
                            end

                            % check reverse direction
                            if model.ub(k) < 0

                                % find reactant
                                reactant_found = false;
                                for m = find(full(model.S(:,k)) > 0)'
                                    if any(m == nad)
                                        reactant_found = true;
                                        break;
                                    end
                                end

                                % find product
                                product_found = false;
                                for m = find(full(model.S(:,k)) < 0)'
                                    if any(m == nadh)
                                        product_found = true;
                                        break;
                                    end
                                end

                                % add
                                if reactant_found && product_found && ~any(k == rxnlist)
                                    rxnlist(end+1) = k;
                                end
                            end
                        end

                    % all nadh --> nad+ reactions
                    case 4

                        % find nad
                        nad = [];
                        for k = 1:length(model.mets)
                            if length(model.mets{k}) == 6
                                if strcmp(model.mets{k}(1:4),'nad[')
                                    nad(end+1) = k;
                                end
                            end
                        end

                        % find nadh
                        nadh = [];
                        for k = 1:length(model.mets)
                            if length(model.mets{k}) == 7
                                if strcmp(model.mets{k}(1:5),'nadh[')
                                    nadh(end+1) = k;
                                end
                            end
                        end

                        % find reactions
                        for k = 1:length(model.rxns)

                            % check forward direction
                            if model.ub(k) > 0

                                % find reactant
                                reactant_found = false;
                                for m = find(full(model.S(:,k)) < 0)'
                                    if any(m == nadh)
                                        reactant_found = true;
                                        break;
                                    end
                                end

                                % find product
                                product_found = false;
                                for m = find(full(model.S(:,k)) > 0)'
                                    if any(m == nad)
                                        product_found = true;
                                        break;
                                    end
                                end

                                % add
                                if reactant_found && product_found && ~any(k == rxnlist)
                                    rxnlist(end+1) = k;
                                end
                            end

                            % check reverse direction
                            if model.ub(k) < 0

                                % find reactant
                                reactant_found = false;
                                for m = find(full(model.S(:,k)) > 0)'
                                    if any(m == nadh)
                                        reactant_found = true;
                                        break;
                                    end
                                end

                                % find product
                                product_found = false;
                                for m = find(full(model.S(:,k)) < 0)'
                                    if any(m == nad)
                                        product_found = true;
                                        break;
                                    end
                                end

                                % add
                                if reactant_found && product_found && ~any(k == rxnlist)
                                    rxnlist(end+1) = k;
                                end
                            end
                        end


                    % all exchange reactions
                    case 5

                        % find reactions
                        for k = 1:length(model.rxns)
                            if any(strcmp(exchange_rxn,model.rxns{k})) && ~any(k == rxnlist)
                                rxnlist(end+1) = k;
                            end
                        end 
                end
            end

            % implement objective value constraint
            model.S(end+1,:) = zeros(1,length(model.rxns));
            model.S(end,find(model.c)) = 1;
            model.b(end+1) = fva_select_fraction * result.objval;

            % clear model.c
            model.c = zeros(size(model.c));

            % calculate FVA fluxes
            flux_min = nan(length(model.rxns),1);
            flux_max = nan(length(model.rxns),1);
            for j = rxnlist

                % set model.c
                model.c(j) = 1;

                % min value
                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = [repmat('=',1,length(model.mets)),'>'];
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'min';

                params.outputflag = 0;
                result = gurobi(model,params);
                try
                    flux_min(j) = result.objval;
                catch
                    flux_min(j) = nan;
                end

                % max value
                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = [repmat('=',1,length(model.mets)),'>'];
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                try
                    flux_max(j) = result.objval;
                catch
                    flux_max(j) = nan;
                end

                % clear model.c
                model.c = zeros(size(model.c));
            end

            % output results
            f_fva = fopen(sprintf('results/%s/%s/fva/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
            fprintf(f_fva,'REACTION\tMIN [mmol/gDW/hr]\tMAX [mmol/gDW/hr]\tLB [mmol/gDW/hr]\tUB [mmol/gDW/hr]\tNAME\tFORMULA\tSUBSYSTEM\tGR RULE\tEC NUMBER\tKEGG\tNOTES\tREFERENCES\n');
            for j = rxnlist
                fprintf(f_fva,'%s\t%0.9f\t%0.9f\t%0.9f\t%0.9f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',model.rxns{j},flux_min(j),flux_max(j),model.lb(j),model.ub(j),model.rxnNames{j},model.rxnFormulas{j},model.subSystems{j},model.grRules{j},model.rxnECNumbers{j},model.rxnKEGGID{j},model.rxnNotes{j},model.rxnReferences{j});
            end
            fclose(f_fva);

        % uniform sampling
        elseif type_of_analysis == 5

            % not yet implemented
            error('Uniform sampling not yet implemented!')

        % media sensitivity - critical value
        elseif type_of_analysis == 6

            % ensure sum(model.c) = 1
            if sum(model.c) ~= 1
                error('ERROR - General - Media Sensitivity - Critical Value - Can only perform media sensitivity with one objective function reaction')
            end
            
            % ensure only one maximization objective
            if length(objective_function_reaction) ~= 1
                error('ERROR - General - Type of Analysis - Can only perform media sensitivity with maximizing one objective function reaction');
            elseif ~strcmp(objective_function_direction{1},'MAX')
                error('ERROR - General - Type of Analysis - Can only perform media sensitivity with maximizing one objective function reaction');
            end

            % check metabolite 1 name
            if ~any(strcmp(model.mets,sprintf('%s[e]',media_sensitivity_critical_metabolite1)))
                error('ERROR - General - Media Sensitivity Options - Critical Value - Metabolite - Metabolite 1 not found in extracellular compartment')

            % check if exchange reaction exists
            elseif ~any(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_critical_metabolite1)))
                error('ERROR - General - Media Sensitivity Options - Critical Value - Metabolite - Metabolite 1 not involved in exchange reaction')

            % ensure that exchange reaction is in exchange_rxn
            elseif ~any(strcmp(exchange_rxn,sprintf('EX_%s[e]',media_sensitivity_critical_metabolite1)))
                error('ERROR - General - Media Sensitivity Options - Critical Value - Metabolite - Exchange reaction for metabolite 1 not found in media composition')
            end

            % find minimum exchange value while maintaining objective value
            original_objective = find(model.c);
            original_objective_value = result.objval;
            model.lb(original_objective) = media_sensitivity_critical_fraction*original_objective_value;
            model.c = zeros(size(model.c));
            model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_critical_metabolite1))) = 1;

            model.A = model.S;
            model.rhs = model.b;
            model.obj = model.c;
            model.sense = repmat('=',1,length(model.mets));
            model.vtype = repmat('C',1,length(model.rxns));
            model.varnames = model.rxns;
            model.modelsense = 'max';

            params.outputflag = 0;
            result = gurobi(model,params);
            critical_val = result.objval;

            % if critical value is negative
            if critical_val < 0

                % save negative of value
                fprintf(f_media,'%s\t%s\t%f\n',samples_all{i},samples_all_source{i},-critical_val);

            % if critical value is positive or zero
            else

                % save nan
                fprintf(f_media,'%s\t%s\tnan\n',samples_all{i},samples_all_source{i});
            end

        % media sensitivity - spectrum
        elseif type_of_analysis == 7

            % ensure sum(model.c) = 1
            if sum(model.c) ~= 1
                error('ERROR - General - Media Sensitivity - Can only perform media sensitivity with one objective function reaction')

            % ensure only one maximization objective
            if length(objective_function_reaction) ~= 1
                error('ERROR - General - Type of Analysis - Can only perform media sensitivity with maximizing one objective function reaction');
            elseif ~strcmp(objective_function_direction{1},'MAX')
                error('ERROR - General - Type of Analysis - Can only perform media sensitivity with maximizing one objective function reaction');
            end

            % check metabolite 1 name
            if ~any(strcmp(model.mets,sprintf('%s[e]',media_sensitivity_spectrum_metabolite1)))
                error('ERROR - General - Media Sensitivity Options - Metabolite - Metabolite 1 not found in extracellular compartment')

            % check if exchange reaction exists
            elseif ~any(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)))
                error('ERROR - General - Media Sensitivity Options - Metabolite - Metabolite 1 not involved in exchange reaction')

            % ensure that exchange reaction is in exchange_rxn
            elseif ~any(strcmp(exchange_rxn,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)))
                error('ERROR - General - Media Sensitivity Options - Metabolite - Exchange reaction for metabolite 1 not found in media composition')
            end  

            if ~isnan(media_sensitivity_spectrum_metabolite2)

                % check metabolite 2 name
                if ~any(strcmp(model.mets,sprintf('%s[e]',media_sensitivity_spectrum_metabolite2)))
                    error('ERROR - General - Media Sensitivity Options - Metabolite - Metabolite 2 not found in extracellular compartment')

                % check if exchange reaction exists
                elseif ~any(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)))
                    error('ERROR - General - Media Sensitivity Options - Metabolite - Metabolite 2 not involved in exchange reaction')

                % ensure that exchange reaction is in exchange_rxn
                elseif ~any(strcmp(exchange_rxn,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)))
                    error('ERROR - General - Media Sensitivity Options - Metabolite - Exchange reaction for metabolite 2 not found in media composition')
                end
            end

            % if 1 metabolite
            if isnan(media_sensitivity_spectrum_metabolite2)

                % find minimum exchange value while maintaining objective value
                original_objective = find(model.c);
                original_objective_value = result.objval;
                original_objective_lb = model.lb(original_objective);
                model.lb(original_objective) = media_sensitivity_spectrum_fraction*original_objective_value;
                model.c = zeros(size(model.c));
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 1;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                start_val = result.objval;

                % if critical value is negative
                if start_val < 0

                    % return original objective function
                    model.c = zeros(size(model.c));
                    model.c(original_objective) = 1;
                    model.lb(original_objective) = original_objective_lb;

                    % initiailize output file
                    f_media = fopen(sprintf('results/%s/%s/media_sensitivity/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
                    fprintf(f_media,'%s UPTAKE LB [mmol/gDW/hr]\tOBJVAL\n',sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1));

                    % progressively increase flux lb
                    uptake = [];
                    objval = [];
                    for j = 1:media_sensitivity_spectrum_points

                        % set lower bound
                        model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = start_val*(media_sensitivity_points-j)/(media_sensitivity_points-1);

                        % find objective value
                        model.A = model.S;
                        model.rhs = model.b;
                        model.obj = model.c;
                        model.sense = repmat('=',1,length(model.mets));
                        model.vtype = repmat('C',1,length(model.rxns));
                        model.varnames = model.rxns;
                        model.modelsense = 'max';

                        params.outputflag = 0;
                        result = gurobi(model,params);

                        % output objective value
                        fprintf(f_media,'%f\t%f\n',model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))),result.objval);
                        uptake(end+1) = -model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)));
                        objval(end+1) = result.objval;
                    end

                    % close output file
                    fclose(f_media);

                    % create image
                    figure; hold on;
                    plot(uptake,objval,'k-');
                    plot([0,max(uptake)],[max(objval),max(objval)],'k--');
                    hold off;
                    xlabel(sprintf('%s Uptake Rate [mmol/gDW/hr]',media_sensitivity_spectrum_metabolite1),'Interpreter','none');
                    ylabel('Objective Value');
                    title(sprintf('Critical %s Uptake Rate = %g mmol/gDW/hr',media_sensitivity_spectrum_metabolite1,max(uptake)),'Interpreter','none');
                    saveas(gcf,sprintf('results/%s/%s/media_sensitivity/%s.png',output_folder,input_filename,samples_all{i}));        

                % if critical value is positive or zero
                else

                    % create output file
                    f_media = fopen(sprintf('results/%s/%s/media_sensitivity/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
                    fprintf(f_media,'%s UPTAKE LB [mmol/gDW/hr]\tOBJVAL\n',sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1));    
                    fprintf(f_media,'nan\t%f\n',original_objective_value);
                    fclose(f_media);          
                end

            % if 2 metabolites
            else

                % find minimum exchanges value while maintaining objective value
                original_objective = find(model.c);
                original_objective_value = result.objval;
                original_objective_lb = model.lb(original_objective);
                model.lb(original_objective) = original_objective_value;
                model.c = zeros(size(model.c));

                % % 1-1 ratio
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 1;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = 1;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                start_val1 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)));
                start_val2 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));

                % % 10-1 ratio
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 10;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = 1;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) < start_val1
                    start_val1 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)));
                end
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) < start_val2
                    start_val2 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));
                end

                % % 1-10 ratio
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 1;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = 10;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) < start_val1
                    start_val1 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)));
                end
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) < start_val2
                    start_val2 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));
                end

                % if both critical values are negative
                if (start_val1 < 0) && (start_val2 < 0)

                    % result original objective function
                    model.c = zeros(size(model.c));
                    model.c(original_objective) = 1;
                    model.lb(original_objective) = original_objective_lb;

                    % initiailize output file
                    f_media = fopen(sprintf('results/%s/%s/media_sensitivity/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
                    fprintf(f_media,'%s UPTAKE LB [mmol/gDW/hr]\t%s UPTAKE LB [mmol/gDW/hr]\tOBJVAL\n',sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1),sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2));

                    % progressively increase flux lb
                    uptake1 = [];
                    uptake2 = [];
                    objval = [];
                    for j = 1:media_sensitivity_points
                        for k = 1:media_sensitivity_points

                            % set lower bound
                            model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = start_val1*(media_sensitivity_points-j)/(media_sensitivity_points-1);
                            model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = start_val2*(media_sensitivity_points-k)/(media_sensitivity_points-1);

                            % find objective value
                            model.A = model.S;
                            model.rhs = model.b;
                            model.obj = model.c;
                            model.sense = repmat('=',1,length(model.mets));
                            model.vtype = repmat('C',1,length(model.rxns));
                            model.varnames = model.rxns;
                            model.modelsense = 'max';

                            params.outputflag = 0;
                            result = gurobi(model,params);

                            % output objective value
                            fprintf(f_media,'%f\t%f\t%f\n',model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))),model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))),result.objval);
                            uptake1(end+1) = -model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)));
                            uptake2(end+1) = -model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));
                            objval(end+1) = result.objval;
                        end
                    end            

                    % close output file
                    fclose(f_media);

                    % create image
                    Z = reshape(objval,media_sensitivity_points,media_sensitivity_points);
                    uptake2 = uptake2;
                    figure;
                    imagesc(uptake1,uptake2,Z);
                    xlabel(sprintf('%s Uptake Rate [mmol/gDW/hr]',media_sensitivity_spectrum_metabolite1),'Interpreter','none');
                    ylabel(sprintf('%s Uptake Rate [mmol/gDW/hr]',media_sensitivity_spectrum_metabolite2),'Interpreter','none');
                    c = colorbar;
                    c.Label.String = 'Objective Value';
                    c.Label.FontSize = 12;
                    set(gca,'YDir','normal');
                    saveas(gcf,sprintf('results/%s/%s/media_sensitivity/%s.png',output_folder,input_filename,samples_all{i}));

                % if one critical value is not negative
                else

                    % create output file
                    f_media = fopen(sprintf('results/%s/%s/media_sensitivity/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
                    fprintf(f_media,'%s UPTAKE LB [mmol/gDW/hr]\t%s UPTAKE LB [mmol/gDW/hr]\tOBJVAL\n',sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1),sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2));
                    fprintf(f_media,'nan\tnan\t%f\n',original_objective_value);
                    fclose(f_media);
                end
            end
        end
    end
        
    % close objective value file
    %fclose(f_obj);
    
    % close critical media sensitivity file
    if type_of_analysis == 6
        fclose(f_media);
    end

end