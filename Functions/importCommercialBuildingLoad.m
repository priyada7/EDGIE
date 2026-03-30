%% NEW FUNCTION: Import Commercial Building Load Data
function [pWorkBase] = importCommercialBuildingLoad(commercial_data_dir,stateName, countyName,weatherTime)
    % Import commercial building electricity consumption data for a specific county
    % and extract data for winter and summer periods

    % Construct directory path for the state
    county_dir = fullfile(commercial_data_dir, stateName);

    % Build the expected filename
    expected_file = [countyName '.csv'];
    county_file = fullfile(county_dir, expected_file);

    % If the exact file doesn't exist, try to find a match
    if exist(county_file, 'file') ~= 2
        files = dir(fullfile(county_dir, '*.csv'));
        match_idx = find(contains({files.name}, countyName, 'IgnoreCase', true), 1);

        if ~isempty(match_idx)
            county_file = fullfile(county_dir, files(match_idx).name);
        else
            fprintf('Warning: Commercial building data not found for %s, %s. Using zeros.\n', countyName, stateName);
            return;
        end
    end


        % Load CSV
       % fprintf('Loading commercial data from: %s\n', county_file);
        commercial_data = readtable(county_file);

        % Validate timestamp column
        if ~any(strcmp(commercial_data.Properties.VariableNames, 'timestamp'))
            fprintf('Warning: No timestamp column found in %s. Using zeros.\n', county_file);
            return;
        end

        % Ensure timestamp is datetime
        if ~isdatetime(commercial_data.timestamp)
            try
                commercial_data.timestamp = datetime(commercial_data.timestamp);
            catch
                fprintf('Warning: Could not convert timestamps in %s. Using zeros.\n', county_file);
                return;
            end
        end

        % Identify building load columns (everything except timestamp)
        building_columns = setdiff(commercial_data.Properties.VariableNames, {'timestamp'});

        if isempty(building_columns)
            fprintf('Warning: No building data columns found in %s. Using zeros.\n', county_file);
            return;
        end

        % Sum loads across all building columns
        total_commercial_load = sum(table2array(commercial_data(:, building_columns)), 2, 'omitnan');

        %% Align years with simulation
      
 

Power = total_commercial_load;                        % individual power profiles
CommercialPower= timetable(commercial_data.timestamp,Power);
CommercialPower = retime(CommercialPower,weatherTime);
pWorkBase=CommercialPower.Var1;



      


    

end