function [cooling_val_95F, cooling_val_82F, heating_val_5C, heating_val_45C] = select_best_hvac_Cooling(x_target, y_target)
    
   % Define available files
    file_list = {'2ton_A_cooling.csv', '2ton_B_cooling.csv', ...
        '3ton_C_cooling.csv', '3ton_A_cooling.csv', '3ton_B_cooling.csv', ...
        '4ton_C_cooling.csv', '4ton_A_cooling.csv', '4ton_B_cooling.csv', ...
        '5ton_C_cooling.csv', '5ton_A_cooling.csv', '5ton_B_cooling.csv'};

   
 
    % Initialize best cooling unit selection
    best_curve_error = inf;
    best_tonnage = '';
    best_manufacturer = '';
    best_cooling_data = [];
    
    % --- COOLING SELECTION ---
    for i = 1:length(file_list)
        filename = file_list{i};
        tokens = regexp(filename, '([2-5]ton)_([a-zA-Z]+)_cooling\.csv', 'tokens');
        
        if isempty(tokens)
            continue;
        end
        
        tonnage_str = tokens{1}{1};
        manufacturer_str = tokens{1}{2};
        
        try
            tbl = readtable(filename);
            x = f2c(tbl{:, 1});
            y = btu2wh(tbl{:, 2});
            
            % Linear fit
            coeffs = polyfit(x, y, 1);
            y_pred = polyval(coeffs, x_target);
            error = abs(y_pred - y_target);
            
            if error < best_curve_error
                best_curve_error = error;
                best_cooling_file = filename;
                best_tonnage = tonnage_str;
                best_manufacturer = manufacturer_str;
                best_cooling_data = coeffs;
            end
        catch
            fprintf('Error loading %s\n', filename);
        end
    end
    
    if isempty(best_cooling_data)
        error('No valid cooling data found.');
    end
    
    % Compute cooling values
    cooling_val_82F = polyval(best_cooling_data, f2c(82));
    cooling_val_95F = polyval(best_cooling_data, f2c(95));
      %   fprintf('Best Cooling Unit: %s\n', best_cooling_file);
    % --- HEATING SELECTION ---
    heating_filename = sprintf('%s_%s_heating.csv', best_tonnage, best_manufacturer);
    
    try
        tbl = readtable(heating_filename);
        x = f2c(tbl{:, 1});
        y = btu2wh(tbl{:, 2});
        
        heating_val_5C = interp1(x, y, f2c(5), 'linear', 'extrap');
        heating_val_45C = interp1(x, y, f2c(45), 'linear', 'extrap');
        %    fprintf('Best Heating Unit: %s\n', heating_filename);
    catch
        error('Error loading heating data for %s.', heating_filename);
    end
 
 
end
