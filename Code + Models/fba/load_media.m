
function [model,concentration_mets,concentration_values,exchange_rxn] = load_media(model,media_choice,media_constraint_metabolite,media_constraint_uptake,concentration_mets,concentration_values,module_5fu,module_cis,module_cpa,module_dox)

    % get list of available media files
    media_dir = dir('../data/media/output/*.csv');
    available_media = {};
    for i = 1:length(media_dir)
        available_media{end+1} = media_dir(i).name;
    end

    % check if media files exists
    for i = 1:length(media_choice)
        if ~any(strcmp(available_media,sprintf('%s.csv',media_choice{i})))
            error('ERROR - Can''t find media file %s.csv',media_choice{i})
        end
    end

    % initialize exchange reactions and flux values
    exchange_rxn = {};
    exchange_value = [];

    % iterate over media files
    for i = 1:length(media_choice)

        % load media file
        f = fopen(sprintf('../data/media/output/%s.csv',media_choice{i}),'r');
        data = textscan(f,'%s %s %f %f','Delimiter',',','headerLines',1);
        fclose(f);

        % add to exchange reactions
        for j = 1:length(data{2})
            if ~any(strcmp(exchange_rxn,data{2}{j}))
                exchange_rxn{end+1} = data{2}{j};
            end
        end
        
        % add to concentration list
        for j = 1:length(data{1})
            if ~any(strcmp(concentration_mets,data{1}{j}))
                concentration_mets{end+1} = data{1}{j};
                concentration_values(end+1) = data{3}(j);
            else
                concentration_values(strcmp(concentration_mets,data{1}{j})) = concentration_values(strcmp(concentration_mets,data{1}{j})) + data{3}(j);
            end
        end
        
    end

    % add exchange reactions for oxygen, carbon dioxide, water, protons, hydroxyl
    exchange_rxn = [exchange_rxn,{'EX_o2[e]','EX_co2[e]','EX_h2o[e]','EX_h[e]','EX_oh1[e]'}];
    
    % add exchange reactions for drug modules
    if module_5fu == 1
        exchange_rxn = [exchange_rxn,{'EX_5fu[e]'}];
    end
    if module_cis == 1
        exchange_rxn = [exchange_rxn,{'EX_cis[e]'}];
    end
    if module_cpa == 1
        exchange_rxn = [exchange_rxn,{'EX_cpa[e]'}];
    end
    if module_dox == 1
        exchange_rxn = [exchange_rxn,{'EX_doxQ[e]'}];
    end
    
    % parse specific media constraints
    media_constraint_exchange = {};
    for i = 1:length(media_constraint_metabolite)
        if ~any(strcmp(model.mets,sprintf('%s[e]',media_constraint_metabolite{i})))
            error('ERROR - Media - Specific Media Constraints - Extracellular metabolite %s not found',media_constraint_metabolite{i})
        elseif ~any(strcmp(model.rxns,sprintf('EX_%s[e]',media_constraint_metabolite{i})))
            error('ERROR - Media - Specific Media Constraints - Exchange reaction for metabolite %s not found',media_constraint_metabolite{i})
        elseif ~any(strcmp(exchange_rxn,sprintf('EX_%s[e]',media_constraint_metabolite{i})))
            error('ERROR - Media - Specific Media Constraints - Metabolite %s not found in selected medium makeup',media_constraint_metabolite{i})
        else
            media_constraint_exchange{end+1} = sprintf('EX_%s[e]',media_constraint_metabolite{i});
        end
    end
    
    % restrict input exchange reactions to those in media
    for i = 1:length(model.rxns)
        if length(model.rxns{i}) >= 3
            if strcmp(model.rxns{i}(1:3),'EX_')
                if ~any(strcmp(exchange_rxn,model.rxns{i}))
                    model.lb(i) = 0;
                else
                    if any(strcmp(media_constraint_exchange,model.rxns{i}))
                        model.lb(i) = -media_constraint_uptake(strcmp(media_constraint_exchange,model.rxns{i}));
                    end
                end
            end
        end
    end
end