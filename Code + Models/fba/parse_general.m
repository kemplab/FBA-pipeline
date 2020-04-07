
function [type_of_analysis,pfba_fraction,fva_all_fraction,fva_select_fraction,fva_select_groups,fva_select_reactions,sampler_fraction,sampler_samples,sampler_steps,media_sensitivity_critical_fraction,media_sensitivity_critical_metabolite1,media_sensitivity_spectrum_fraction,media_sensitivity_spectrum_metabolite1,media_sensitivity_spectrum_metabolite2,media_sensitivity_spectrum_points] = parse_general(input_file)

    % load data
    [~,~,raw] = xlsread(input_file,'General');

    % type of analysis
    for row = 2:8
        value = raw(row,1);
        value = value{1};
        if ~isnan(value)
            if ischar(value)
                if strcmpi(value,'X')
                    if exist('type_of_analysis','var')
                        error('ERROR - General - Type of Analysis - Can only choose one type of analysis')
                    else
                        type_of_analysis = row - 1;
                    end
                else
                    error('ERROR - General - Type of Analysis - Values must either be X or blank')
                end
            else
                error('ERROR - General - Type of Analysis - Values must either be X or blank')
            end
        end
    end
    if ~exist('type_of_analysis','var')
        error('ERROR - General - Type of Analysis - Must select one type of analysis')
    end
    
    % options - pfba
    if type_of_analysis == 2
    
        % fraction of maximum objective value
        value = raw(11,1);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)
                if (value >= 0) && (value <= 1)
                    pfba_fraction = value;
                else
                    error('ERROR - General - Options - pFBA - Fraction of maximum objective value must either be >= 0 and <= 1')
                end
            else
                error('ERROR - General - Options - pFBA - Fraction of maximum objective value must either be >= 0 and <= 1')
            end
        else
            error('ERROR - General - Options - pFBA - Must enter fraction of maximum objective value')
        end
        if ~exist('pfba_fraction','var')
            error('ERROR - General - Options - pFBA - Must enter fraction of maximum objective value')
        end 
    
    else
        pfba_fraction = NaN; 
    end
    
    % options - fva - all reactions
    if type_of_analysis == 3
    
        % fraction of maximum objective value
        value = raw(14,1);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)
                if (value >= 0) && (value <= 1)
                    fva_all_fraction = value;
                else
                    error('ERROR - General - Options - FVA - All reactions - Fraction of maximum objective value must either be >= 0 and <= 1')
                end
            else
                error('ERROR - General - Options - FVA - All reactions - Fraction of maximum objective value must either be >= 0 and <= 1')
            end
        else
            error('ERROR - General - Options - FVA - All reactions - Must enter fraction of maximum objective value')
        end
        if ~exist('fva_all_fraction','var')
            error('ERROR - General - Options - FVA - All reactions - Must enter fraction of maximum objective value')
        end 
        
    else    
        fva_all_fraction = NaN;      
    end
    
    % options - fva - select reactions
    if type_of_analysis == 4
        
        % fraction of maximum objective value
        value = raw(17,1);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)
                if (value >= 0) && (value <= 1)
                    fva_select_fraction = value;
                else
                    error('ERROR - General - Options - FVA - Select reactions - Fraction of maximum objective value must either be >= 0 and <= 1')
                end
            else
                error('ERROR - General - Options - FVA - Select reactions - Fraction of maximum objective value must either be >= 0 and <= 1')
            end
        else
            error('ERROR - General - Options - FVA - Select reactions - Must enter fraction of maximum objective value')
        end
        if ~exist('fva_select_fraction','var')
            error('ERROR - General - Options - FVA - Select reactions - Must enter fraction of maximum objective value')
        end
        
        % reaction groups
        fva_select_groups = [];
        for row = 19:23
            value = raw(row,1);
            value = value{1};
            if ~isnan(value)
                if ischar(value)
                    if strcmpi(value,'X')
                        fva_select_groups(end+1) = row - 18;
                    else
                        error('ERROR - General - Options - FVA - Select groups - Values must either be X or blank')
                    end
                else
                    error('ERROR - General - Options - FVA - Select groups - Values must either be X or blank')
                end
            end
        end

        % reaction list
        fva_select_reactions = {};
        for row = 25:size(raw,1)
            value = raw(row,2);
            value = value{1};
            if ~isnan(value)
                if ischar(value)
                    fva_select_reactions{end+1} = value;
                else
                    error('ERROR - General - Options - FVA - Select reactions - Reaction name must be a string')
                end
            end
        end
        if (length(fva_select_groups) == 0) && (length(fva_select_reactions) == 0)
            error('ERROR - General - Options - FVA - Select reactions - Must provide at least one reaction group or name')
        end
        
    else
        fva_select_fraction = NaN;
        fva_select_groups = [];
        fva_select_reactions = {};
    end
    
    % options - uniform sampling
    if type_of_analysis == 5
        
        % fraction of maximum objective value
        value = raw(11,5);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)
                if (value >= 0) && (value <= 1)
                    sampler_fraction = value;
                else
                    error('ERROR - General - Options - Uniform sampling - Fraction of maximum objective value must either be >= 0 and <= 1')
                end
            else
                error('ERROR - General - Options - Uniform sampling - Fraction of maximum objective value must either be >= 0 and <= 1')
            end
        else
            error('ERROR - General - Options - Uniform sampling - Must enter fraction of maximum objective value')
        end
        if ~exist('sampler_fraction','var')
            error('ERROR - General - Options - Uniform sampling - Must enter fraction of maximum objective value')
        end
        
        % number of samples
        value = raw(12,5);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)
                if mod(value,1) == 0
                    if value >= 1
                        sampler_samples = value;
                    else
                        error('ERROR - General - Options - Uniform sampling - Number of Samples - Value must be >= 1')
                    end
                else
                    error('ERROR - General - Options - Uniform sampling - Number of Samples - Value must be an integer')
                end
            else
                error('ERROR - General - Options - Uniform sampling - Number of Samples - Value must be a number')
            end
        end
        if ~exist('sampler_samples','var')
            error('ERROR - General - Options - Uniform sampling - Number of Samples - Must enter value')
        end

        % number of steps per sample
        value = raw(13,5);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)
                if mod(value,1) == 0
                    if value >= 1
                        sampler_steps = value;
                    else
                        error('ERROR - General - Options - Uniform sampling - Number of Steps per Sample - Value must be >= 1')
                    end
                else
                    error('ERROR - General - Options - Uniform sampling - Number of Steps per Sample - Value must be an integer')
                end
            else
                error('ERROR - General - Options - Uniform sampling - Number of Steps per Sample - Value must be a number')
            end
        end
        if ~exist('sampler_steps','var')
            error('ERROR - General - Options - Uniform sampling - Number of Steps per Sample - Must enter value')
        end

    else
        sampler_fraction = NaN;
        sampler_samples = NaN;
        sampler_steps = NaN;
    end
    
    % options - media sensitivity - critical value
    if type_of_analysis == 6
        
        % fraction of maximum objective value
        value = raw(16,5);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)
                if (value >= 0) && (value <= 1)
                    media_sensitivity_critical_fraction = value;
                else
                    error('ERROR - General - Options - Media sensitivity - Critical value - Fraction of maximum objective value must either be >= 0 and <= 1')
                end
            else
                error('ERROR - General - Options - Media sensitivity - Critical value - Fraction of maximum objective value must either be >= 0 and <= 1')
            end
        else
            error('ERROR - General - Options - Media sensitivity - Critical values - Must enter fraction of maximum objective value')
        end
        if ~exist('media_sensitivity_critical_fraction','var')
            error('ERROR - General - Options - Media sensitivity - Critical value - Must enter fraction of maximum objective value')
        end
        
        % metabolite 1
        value = raw(17,5);
        value = value{1};
        if ~isnan(value)
            if ischar(value)
                media_sensitivity_critical_metabolite1 = value;
            else
                error('ERROR - General - Options - Media sensitivity - Critical value - Metabolite - Metabolite 1 name must be a string')
            end 
        end
        if ~exist('media_sensitivity_critical_metabolite1','var')
            error('ERROR - General - Options - Media sensitivity - Critical value - Metabolite - Must enter metabolite 1 name')
        end
    
    else
        media_sensitivity_critical_fraction = NaN;
        media_sensitivity_critical_metabolite1 = NaN;    
    end
    
    % options - media sensitivity - spectrum
    if type_of_analysis == 7
        
        % fraction of maximum objective value
        value = raw(20,5);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)
                if (value >= 0) && (value <= 1)
                    media_sensitivity_spectrum_fraction = value;
                else
                    error('ERROR - General - Options - Media sensitivity - Spectrum - Fraction of maximum objective value must either be >= 0 and <= 1')
                end
            else
                error('ERROR - General - Options - Media sensitivity - Spectrum - Fraction of maximum objective value must either be >= 0 and <= 1')
            end
        else
            error('ERROR - General - Options - Media sensitivity - Spectrum - Must enter fraction of maximum objective value')
        end
        if ~exist('media_sensitivity_spectrum_fraction','var')
            error('ERROR - General - Options - Media sensitivity - Spectrum - Must enter fraction of maximum objective value')
        end
        
        % metabolite 1
        value = raw(21,5);
        value = value{1};
        if ~isnan(value)
            if ischar(value)
                media_sensitivity_spectrum_metabolite1 = value;
            else
                error('ERROR - General - Options - Media sensitivity - Spectrum - Metabolite - Metabolite 1 name must be a string')
            end 
        end
        if ~exist('media_sensitivity_spectrum_metabolite1','var')
            error('ERROR - General - Options - Media sensitivity - Spectrum - Metabolite - Must enter metabolite 1 name')
        end
        
        % metabolite 2
        value = raw(22,5);
        value = value{1};
        if ~isnan(value)
            if ischar(value)
                media_sensitivity_spectrum_metabolite2 = value;
            else
                error('ERROR - General - Options - Media sensitivity - Spectrum - Metabolite - Metabolite 2 name must be a string')
            end 
        end
        if ~exist('media_sensitivity_spectrum_metabolite2','var')
            media_sensitivity_spectrum_metabolite2 = NaN;
        end
        
        % number of concentration points
        value = raw(23,5);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)
                if mod(value,1) == 0
                    if value >= 1
                        media_sensitivity_spectrum_points = value;
                    else
                        error('ERROR - General - Options - Media sensitivity - Spectrum - Number of Concentration Points - Value must be >= 1')
                    end
                else
                    error('ERROR - General - Options - Media sensitivity - Spectrum - Number of Concentration Points - Value must be an integer')
                end
            else
                error('ERROR - General - Options - Media sensitivity - Spectrum - Number of Concentration Points - Value must be a number')
            end
        end
        if ~exist('media_sensitivity_spectrum_points','var')
            error('ERROR - General - Options - Media sensitivity - Spectrum - Number of Concentration Points - Must enter value')
        end
    
    else
        media_sensitivity_spectrum_fraction = NaN;
        media_sensitivity_spectrum_metabolite1 = NaN;
        media_sensitivity_spectrum_metabolite2 = NaN;
        media_sensitivity_spectrum_points = NaN;
    end  
end