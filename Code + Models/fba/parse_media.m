
function [media_choice,media_constraint_metabolite,media_constraint_uptake] = parse_media(input_file)

    % load data
    [~,~,raw] = xlsread(input_file,'Media');
    
    % media to add
    media_choice = {};
    for row = 9:size(raw,1)
        value = raw(row,1);
        value = value{1};

        if ~isnan(value)
            if ischar(value)
                media_choice{end+1} = value;
            else
                error('ERROR - Media - Medium Makeup - Value must either be X or blank')
            end
        end
    end
    
    % specific media constraints
    media_constraint_metabolite = {};
    media_constraint_uptake = [];

    for row = 3:size(raw,1)
        value = raw(row,5);
        value = value{1};
        
        % metabolite name
        if ~isnan(value)
            if ischar(value)
                media_constraint_metabolite{end+1} = value;

                % uptake rate
                value = raw(row,6);
                value = value{1};
                if ~isnan(value)
                    if isfloat(value)           
                        if value >= 0
                            media_constraint_uptake(end+1) = value;
                        else
                            error('ERROR - Media - Specific Media Constraints - Maximum uptake rate must be >= 0')
                        end
                    else
                        error('ERROR - Media - Specific Media Constraints - Maximum uptake rate must be a number')
                    end
                else
                    error('ERROR - Media - Specific Media Constraints - Must give maximum uptake rate')
                end
            else
                error('ERROR - Media - Specific Media Constraints - Metabolite name must be a string')
            end
        end
    end
end