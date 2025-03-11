function [cooling_val_95F, cooling_val_82F, heating_val_5C, heating_val_45C] = select_best_hvac_Cooling(x_target, y_target)
    
   % Define available files
    file_list = {'2ton_A_cooling.csv', '2ton_B_cooling.csv', ...
        '3ton_C_cooling.csv', '3ton_A_cooling.csv', '3ton_B_cooling.csv', ...
        '4ton_C_cooling.csv', '4ton_A_cooling.csv', '4ton_B_cooling.csv', ...
        '5ton_C_cooling.csv', '5ton_A_cooling.csv', '5ton_B_cooling.csv',...
        '6ton_C_cooling.csv', '6ton_A_cooling.csv', '6ton_B_cooling.csv',...
        '7ton_C_cooling.csv', '7ton_A_cooling.csv', '7ton_B_cooling.csv',...
        '8ton_C_cooling.csv', '8ton_A_cooling.csv', '8ton_B_cooling.csv'};
 
%    
%   file_list = {'3ton_C_cooling.csv', '3ton_A_cooling.csv', '3ton_B_cooling.csv'};
    % Initialize best cooling unit selection
    best_curve_error = inf;
    best_tonnage = '';
    best_manufacturer = '';
    best_cooling_data = [];
    
%     figure;
%     hold on;
%     xlabel('Temperature (°C)');
%     ylabel('Cooling Capacity (Wh)');
%     title('Linear Fits of HVAC Cooling Units');
%     grid on;
%     colors = lines(length(file_list));  % Generate different colors for each curve
    
    % --- COOLING SELECTION ---
    for i = 1:length(file_list)
        filename = file_list{i};
        tokens = regexp(filename, '([2-8]ton)_([a-zA-Z]+)_cooling\.csv', 'tokens');
        
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
            
            % Linear fit
%            % coeffs = polyfit(x, y, 1);
%             x_fit = linspace(min(x), max(x), 100);
%             y_fit = polyval(coeffs, x_fit);
%             
%             % Plot each linear fit
%             plot(x_fit, y_fit, 'Color', colors(i, :), 'LineWidth', 1.5, 'DisplayName', filename);
%             
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
       % fprintf('Best Cooling Unit: %s\n', best_cooling_file);
        % Plot target point
   % scatter((x_target), y_target, 100, 'k', 'filled', 'DisplayName', 'Target Point');
    % Plot best-fit line
%     x_best_fit = linspace(f2c(60), f2c(110), 100);  % Extend range for visualization
%     y_best_fit = polyval(best_cooling_data, x_best_fit);
   % plot(x_best_fit, y_best_fit, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Best Fit');
    
%     legend('show');
%     hold off;
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
