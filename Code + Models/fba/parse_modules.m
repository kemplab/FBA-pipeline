
function [catalase_version, module_5fu, module_cis, module_cpa, module_dox] = parse_modules(input_file)

    % load data
    [~,~,raw] = xlsread(input_file,'Modules');
    
    % catalase version
    value = raw(3,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                catalase_version = 1;
            else
                error('ERROR - Modules - Catalase Version - Value must either be X or blank')
            end
        else
            error('ERROR - Modules - Catalase Version - Value must either be X or blank')
        end
    else
        catalase_version = 0;
    end
    
    value = raw(4,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                if catalase_version == 1;
                   error('ERROR - Modules - Catalase Version - Cannot select both options') 
                end
            else
                error('ERROR - Modules - Catalase Version - Value must either be X or blank')
            end
        else
            error('ERROR - Modules - Catalase Version - Value must either be X or blank')
        end
    elseif catalase_version == 0;
        error('ERROR - Modules - Catalase Version - Must select an option')
    end
    
    % 5-fluorouracil module
    value = raw(8,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                module_5fu = 1;
            else
                error('ERROR - Modules - 5-Fluorouracil - Value must either be X or blank')
            end
        else
            error('ERROR - Modules - 5-Fluorouracil - Value must either be X or blank')
        end
    else
        module_5fu = 0;
    end
    
    % cisplatin module
    value = raw(9,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                module_cis = 1;
            else
                error('ERROR - Modules - Cisplatin - Value must either be X or blank')
            end
        else
            error('ERROR - Modules - Cisplatin - Value must either be X or blank')
        end
    else
        module_cis = 0;
    end
    
    % cyclophosphamide module
    value = raw(10,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                module_cpa = 1;
            else
                error('ERROR - Modules - Cyclophosphamide - Value must either be X or blank')
            end
        else
            error('ERROR - Modules - Cyclophosphamide - Value must either be X or blank')
        end
    else
        module_cpa = 0;
    end
    
    % doxorubicin module
    value = raw(11,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                module_dox = 1;
            else
                error('ERROR - Modules - Doxorubicin - Value must either be X or blank')
            end
        else
            error('ERROR - Modules - Doxorubicin - Value must either be X or blank')
        end
    else
        module_dox = 0;
    end
 
end