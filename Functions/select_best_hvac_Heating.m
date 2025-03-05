function [cooling_val_95F, cooling_val_82F, heating_val_5C, heating_val_45C] = select_best_hvac_Heating(x_target, y_target)
    % Define available files
    heating_files = {'2ton_A_heating.csv', '2ton_B_heating.csv', ...
        '3ton_C_heating.csv', '3ton_A_heating.csv', '3ton_B_heating.csv', ...
        '4ton_C_heating.csv', '4ton_A_heating.csv', '4ton_B_heating.csv', ...
        '5ton_C_heating.csv', '5ton_A_heating.csv', '5ton_B_heating.csv'};
    
    cooling_files = {'2ton_A_cooling.csv', '2ton_B_cooling.csv', ...
        '3ton_C_cooling.csv', '3ton_A_cooling.csv', '3ton_B_cooling.csv', ...
        '4ton_C_cooling.csv', '4ton_A_cooling.csv', '4ton_B_cooling.csv', ...
        '5ton_C_cooling.csv', '5ton_A_cooling.csv', '5ton_B_cooling.csv'};

%      heating_files = {'3ton_C_heating.csv', '3ton_A_heating.csv', '3ton_B_heating.csv'};
%     
%     cooling_files = { '3ton_C_cooling.csv', '3ton_A_cooling.csv', '3ton_B_cooling.csv'};

    % Temperature range for heating
    temp_min_heating = f2c(5);
    temp_max_heating = f2c(41);

    % Initialize best heating unit selection
    best_curve_error = inf;
    best_heating_file = '';
    best_heating_data = [];

    % --- HEATING SELECTION ---
    for i = 1:numel(heating_files)
        filename = heating_files{i};
        try
            tbl = readtable(filename);
            if width(tbl) >= 2
                x = f2c(tbl{:, 1});
                y = btu2wh(tbl{:, 2});
                valid_idx = (x >= temp_min_heating) & (x <= temp_max_heating);
                x = x(valid_idx);
                y = y(valid_idx);
                
                % Linear fit
                if length(x) > 1
                    coeffs = polyfit(x, y, 1);
                    y_pred = polyval(coeffs, x_target);
                    error = abs(y_pred - y_target);
                    
                    % Check if this is the best heating unit
                    if error < best_curve_error
                        best_curve_error = error;
                        best_heating_file = filename;
                        best_heating_data = coeffs;
                    end
                end
            end
        catch
            fprintf('Error loading %s\n', filename);
        end
    end

    % Display selected heating unit
    if ~isempty(best_heating_data)
        heating_val_5C = polyval(best_heating_data, f2c(5));
        heating_val_45C = polyval(best_heating_data, f2c(45));
    else
        fprintf('No valid heating data found.\n');
        return;
    end

    % --- COOLING SELECTION ---
    best_cooling_file = strrep(best_heating_file, 'heating', 'cooling');
    try
        tbl = readtable(best_cooling_file);
        if width(tbl) >= 2
            x = f2c(tbl{:, 1});
            y = btu2wh(tbl{:, 2});
            cooling_val_82F = interp1(x, y, f2c(82), 'linear', 'extrap');
            cooling_val_95F = interp1(x, y, f2c(95), 'linear', 'extrap');
        else
            fprintf('No valid cooling data found.\n');
        end
    catch
        fprintf('Error loading %s\n', best_cooling_file);
    end
    
    % Print best file names
%     fprintf('Best Heating File: %s\n', best_heating_file);
%     fprintf('Best Cooling File: %s\n', best_cooling_file);
end
