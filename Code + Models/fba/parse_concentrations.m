
function [recalculate_thermodynamics,default_concentration_lb,default_concentration_ub,concentration_ranges_all,concentration_ranges_sample_folder,concentration_ranges_sample] = parse_concentrations(input_file)

    % load data
    [~,~,raw] = xlsread(input_file,'Concentrations');

    % recalculate thermodynamics
    row = 3;
    value = raw(row,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                recalculate_thermodynamics = 1;
            else
                error('ERROR - Concentrations - Recalculate Thermodynamics - Value must either be X or blank')
            end
        else
            error('ERROR - Concentrations - Recalculate Thermodynamics - Value must either be X or blank')
        end
    else
        recalculate_thermodynamics = 0;
    end

    row = 4;
    value = raw(row,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                if recalculate_thermodynamics == 1;
                   error('ERROR - Concentrations - Recalculate Thermodynamics - Cannot select both Yes and No') 
                end
            else
                error('ERROR - Concentrations - Recalculate Thermodynamics - Value must either be X or blank')
            end
        else
            error('ERROR - Concentrations - Recalculate Thermodynamics - Value must either be X or blank')
        end
    elseif recalculate_thermodynamics == 0;
        error('ERROR - Concentrations - Recalculate Thermodynamics - Must select either Yes and No')
    end
    
    % default concentration ranges - lower bound
    value = raw(8,1);
    value = value{1};
    if recalculate_thermodynamics == 1
        if ~isnan(value)
            if isfloat(value)
                if value >= 0
                    default_concentration_lb = value;
                else
                    error('ERROR - Concentrations - Default Lower Bound - Value must be >= 0')
                end
            else
                error('ERROR - Concentrations - Default Lower Bound - Value must be a number')
            end
        else
            error('ERROR - Concentrations - Default Lower Bound - Must provide value if recalculating thermodynamics')
        end
    else
        default_concentration_lb = NaN;
    end
    
    % default concentration ranges - upper bound
    value = raw(9,1);
    value = value{1};
    if recalculate_thermodynamics == 1
        if ~isnan(value)
            if isfloat(value)
                if value >= 0
                    if value >= default_concentration_lb
                        default_concentration_ub = value;
                    else
                        error('ERROR - Concentrations - Default Upper Bound - Value must be >= default lower bound')
                    end
                else
                    error('ERROR - Concentrations - Default Upper Bound - Value must be >= 0')
                end
            else
                error('ERROR - Concentrations - Default Upper Bound - Value must be a number')
            end
        else
            error('ERROR - Concentrations - Default Upper Bound - Must provide value if recalculating thermodynamics')
        end
    else
        default_concentration_ub = NaN;
    end
    
    % concentration values/ranges applied to all samples
    concentration_ranges_all = {};
    for row = 13:size(raw,1)
        value = raw(row,2);
        value = value{1};

        if ~isnan(value)
            if ischar(value)
                concentration_ranges_all{end+1} = value;
            else
                error('ERROR - Concentrations - Ranges Applied to All Samples - File name must be valid string')
            end
        end
    end
    
    % concentration values/ranges applied to individual samples
    concentration_ranges_sample_folder = {};
    concentration_ranges_sample = {};
    for row = 13:size(raw,1)

        value = raw(row,4);
        value = value{1};
        if ~isnan(value)
            if ischar(value)
            
                % folder name
                concentration_ranges_sample_folder{end+1} = value;
                
                % sample name
                value = raw(row,5);
                value = value{1};
                if ~isnan(value)
                    if ischar(value)
                        concentration_ranges_sample{end+1} = value;
                    else
                        error('ERROR - Concentrations - Ranges Applied to Individual Samples - Sample name must be valid string')
                    end
                else
                    error('ERROR - Concentrations - Ranges Applied to Individual Samples - Must provide sample name')
                end       
            else
                error('ERROR - Concentrations - Ranges Applied to Individual Samples - Folder name must be valid string')
            end
        end
     end
end