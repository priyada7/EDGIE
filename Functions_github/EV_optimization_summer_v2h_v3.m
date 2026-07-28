function [e,p2v,HousingUnitsload,cvx_status]=...
    EV_optimization_summer_v2h_v3(K,p1baseHPSummer,p1basehpwhSummer,a2,e0,eMax,eMin,tau,ec,dt,pdMax,pcMax,atWork,atHome,...
n2,PSummer,pWorkBase,etac,etad,t,tw,th,onRoad,n1,SummerPeak,TodaysHeadroom,FutureHeadroom,SummerPeakComStock,U_units2bldg,A_attached_bldg,A_detached_home,PdetachedSummer,PattachedSummer,n1d,n1a,cars_attached,cars_detached,ev_start_attached,EnergyCommercial,EnergyResidential,normPrice,v2hEligible,A_hv,C_bh,D_bv,scenario,...
atGrocery,groceryEnergy,atErrands,errandEnergy)

% % Complete EV allocation with realistic probability-based distribution
% and building-level sharing for attached homes
%
% Inputs:
%   Standard optimization inputs (K, loads, EV parameters, etc.)
%   frac_detached: fraction of detached homes (0-1)
%   units_per_building: number of units per attached building (default: 20)
%
% Outputs:
%   e, pChem, cvx_status: optimization results
%   A_detached: sparse mapping for detached (n_detached × n_detached_evs)
%   A_attached: sparse mapping for attached (n_buildings × n_attached_evs)
%   EV_to_home: EV assignment array
%   building_to_home: cell array mapping buildings to homes
%   cars_detached, cars_attached, cars_per_building: allocation details

n_detached = n1d;
n_attached = n1a ;



% --------------------------
% Prepare EV Power Limits
% --------------------------
pChemMin = zeros(K, n2);  % Minimum power (discharge limit)
pChemMax = zeros(K, n2);  % Maximum power (charge limit)
w = zeros(K, n2);         % Energy consumption while driving


valDischarge = pcMax ./ etad;  % Max discharge power 
valCharge    = pcMax .* etac;  % Max charge power 


% When EV is at home: allow both charging and discharging
pChemMin = (pChemMin - atHome .* repmat(valDischarge, K, 1));  % Negative = discharge
pChemMax = pChemMax + atHome .* repmat(valCharge, K, 1);     % Positive = charge

%%
 
   
    % Set non-eligible EVs to no-discharge
    pChemMin(:, ~v2hEligible) = 0;



%%
% When EV is at work: only charging allowed (no V2H at work)
pChemMax = pChemMax + atWork .* repmat(valCharge, K, 1);

% When EV is on road: no charging/discharging (set limits to zero)
pChemMin(onRoad) = 0;
pChemMax(onRoad) = 0;

pChemMin(logical(atGrocery)) = 0;   
pChemMax(logical(atGrocery)) = 0;   

pChemMin(logical(atErrands)) = 0;   
pChemMax(logical(atErrands)) = 0;  
for k=1:K
   w(k,onRoad(k,:)) = ec(onRoad(k,:))/dt;
   
   % grocery run energy
    w(k, atGrocery(k,:)==1) = w(k, atGrocery(k,:)==1) + groceryEnergy(atGrocery(k,:)==1) / dt;

    w(k, atErrands(k,:)==1) = w(k, atErrands(k,:)==1) + errandEnergy(atErrands(k,:)==1) / dt;
end

% --------------------------
% --------------------------
pie = 0.15; % Energy price [$/kWh]
pid = 10;   % Demand charge [$/kW]
DeviationPenalty_EV = EnergyResidential / mean(etac); % Penalty for not returning to initial SOC

etad = repmat(etad, K, 1);    
etac = repmat(etac, K, 1);    
eMin = repmat(eMin, K+1, 1);  
eMax = repmat(eMax, K+1, 1);  
a2   = repmat(a2, K, 1);      
tau  = repmat(tau, K, 1);     

% --------------------------
% CVX Optimization 
% --------------------------
fprintf('Setting up CVX optimization for Summer...\n');

cvx_begin quiet
cvx_precision low
    variables e(K+1,n2) pChem(K,n2)

expressions homeLoadSummer workLoadSummer diffCalc;  % just declare holder first

% Total loads
homeLoadSummer = PSummer...% MEL MFRED
            + sum(p1baseHPSummer,2)... % Heat pump
            + sum(p1basehpwhSummer,2) ...% Heat-pump water heater
            + sum(atHome.*max(etad.*pChem, pChem./etac),2); % EV

workLoadSummer = pWorkBase + sum(atWork .* (pChem ./ etac), 2); % ComStock + EVcharging
diffCalc=pChem(2:K,:) - pChem(1:K-1,:);

minimize( ...
    DeviationPenalty_EV*sum(max(0, e0 - e(K+1,:))) ...   % EV deviation penalty
   + (1/K )*dt*EnergyResidential*sum(homeLoadSummer) ...                          % Energy usage cost
  +   (1/K )*dt*EnergyCommercial*sum(workLoadSummer) ...                          % Energy usage cost
    + 1000*pos(FutureHeadroom*max(homeLoadSummer) - TodaysHeadroom*SummerPeak) ...            % Summer peak violation penalty
    + 1000*pos(FutureHeadroom*max(workLoadSummer) - TodaysHeadroom*SummerPeakComStock) ... 
    + (1/(K-1) )*normPrice*norm(diffCalc, 1) ... 
)  


    subject to
    % --------------------------
    % Constraints
    % --------------------------
    e(1,:) == e0;
    0.2*eMax  <= e <= eMax;
    e(2:K+1,:) == a2.*e(1:K,:) + (1-a2).*tau.*(pChem - w);


    pChemMin <= pChem <= pChemMax;

    % DETACHED: Per-home demand (K × n_detached)
    homeDemand_detached = PdetachedSummer+ ...
                         p1baseHPSummer(:,1:n_detached) + ...
                         p1basehpwhSummer(:,1:n_detached);
    
    % ATTACHED: Per-building demand (K × n_buildings)
    % Aggregate all unit demands within each building
    homeDemand_attached_units =PattachedSummer+ ...
                               p1baseHPSummer(:,n_detached+1:n1) + ...
                               p1basehpwhSummer(:,n_detached+1:n1);  % K × n1a
    
    HousingUnitsload =[homeDemand_detached,homeDemand_attached_units];
 
 if scenario == 1
        % ia) building level for attached + detached
        -(D_bv * (pChem .* etad)') <= C_bh * HousingUnitsload';

    elseif scenario == 2
        % ib) per unit for attached, per home for detached — no cross-unit sharing
        -(A_hv * (pChem .* etad)') <= HousingUnitsload';

    elseif scenario == 3
        % ii) pessimistic — V2H only detached, attached cannot discharge
        A_hv_det = A_hv(1:n_detached, :);
        -(A_hv_det * (pChem .* etad)') <= HousingUnitsload(:, 1:n_detached)';
        pChem(:, ev_start_attached:end) >= 0;
    end

cvx_end
e=full(e);
pChem=full(pChem);
p2v = max(etad .* pChem, pChem ./ etac) .* atHome ...   % V2H: at home
    + (pChem ./ etac) .* atWork;      
% --------------------------
% Results Processing
% --------------------------
fprintf('CVX optimization completed with status: %s\n', cvx_status);

% %
% ==========================================
% POST-ANALYSIS VALIDATION PLOTS
% ==========================================
% 
% --------------------------
% 1. DETACHED HOMES: 5 Random Homes
% --------------------------
% ==========================================
% POST-ANALYSIS VALIDATION PLOTS
% ==========================================
% 
% --------------------------
% 1. DETACHED HOMES: 5 Random Homes
% --------------------------
% fprintf('\n=== Generating Detached Homes Validation Plots ===\n');
% rng("shuffle")
% % Select 5 random detached homes
% n_sample_detached = 4;
% detached_sample_idx = randsample(1:n_detached, n_sample_detached, false);
% 
% 
% for plot_idx = 1:n_sample_detached
%     home_id = detached_sample_idx(plot_idx);
%     
%     % Get home load (MEL + HP + HPWH)
%     home_load = PdetachedSummer(:, home_id) + p1baseHPSummer(:, home_id) + p1basehpwhSummer(:, home_id);
%     
%     % Find EVs associated with this home
%     ev_indices = find(EV_to_home(1:ev_start_attached-1) == home_id);
%     
%     % Sum EV discharge power for this home (positive = discharge)
%     if isempty(ev_indices)
%         ev_discharge = zeros(K, 1);
%     else
%         ev_discharge = sum(max(0, -pChem(1:K, ev_indices)).*etad(1:K, ev_indices), 2);  % negative pChem = discharge
%     end
%     
%     % Create subplot
%     ax = subplot(2, 3, plot_idx);
%     
%     time_vec = 1:K;
%     
%     % Plot 1: Home load (bar) and total EV discharge (line overlay)
%     h_bar = bar(time_vec, home_load, 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'none');
%     h_bar.FaceAlpha = 0.6;
%     ylabel(' Load (kW)');
%     set(gca, 'YColor', 'k');
%     hold on
%     plot(time_vec, ev_discharge, 'm-', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3);    
%    
%     xlabel('Time');
%     grid on; grid minor;
%     
%     % Title with home info
%     title(sprintf('Home %d | %.0f sqft | %d EVs', home_id, size_detached(home_id), cars_detached(home_id)), ...
%          'FontWeight', 'bold');
%     
%     set(ax);
% end
% 
% 
% sgtitle('Detached Homes: Load Profile & EV Discharge Validation', 'FontSize', 14, 'FontWeight', 'bold');
% % 
% % % --------------------------
% % % 2. ATTACHED BUILDINGS: 5 Random Buildings
% % % --------------------------
% fprintf('Generating Attached Buildings Validation Plots...\n');
% 
% % Select 5 random attached buildings
% n_sample_buildings = min(4, n_buildings);
% building_sample_idx = randsample(1:n_buildings, n_sample_buildings, false);
% 
% 
% for plot_idx = 1:n_sample_buildings
%     bldg_id = building_sample_idx(plot_idx);
%     
%     % Get all units in this building
%     units_in_bldg = building_to_home{bldg_id};
%     units_in_bldg = units_in_bldg(units_in_bldg >= n_detached+1 & units_in_bldg <= n1);
%     
%     % Sum building load (all units in building)
%     building_load = sum(PattachedWinter(:, units_in_bldg - n_detached) + ...
%                         p1baseHPWinter(:, units_in_bldg) + ...
%                         p1basehpwhWinter(:, units_in_bldg), 2);
%     
%     % Find EVs associated with this building
%     ev_indices = find(EV_to_home(ev_start_attached:end) == bldg_id);
%     ev_indices = ev_indices + ev_start_attached - 1;  % Convert to actual EV indices
%     
%     % Sum EV discharge power for this building
%     if isempty(ev_indices)
%         ev_discharge = zeros(K, 1);
%     else
%         ev_discharge = sum(max(0, -pChem(1:K, ev_indices).*etad(1:K, ev_indices)), 2);  % negative pChem = discharge
% 
%     end
%     
%     % Create subplot
%     ax = subplot(2, 3, plot_idx);
%     
%     time_vec = 1:K;
%     
%     % Plot: Building load and EV discharge
%    
%     h_bar = bar(time_vec, building_load, 'FaceColor', [0.7 0.9 0.7], 'EdgeColor', 'none');
%     h_bar.FaceAlpha = 0.6;
%     ylabel('Building Load (kW)');
%     set(gca, 'YColor', 'k');
%      hold on
%     plot(time_vec, ev_discharge, 'r-', 'LineWidth', 2.5, 'Marker', 's', 'MarkerSize', 4);
%     xlabel('Time Step');
%     grid on; grid minor;
%     
%     % Title with building info
%     n_units = length(units_in_bldg);
%     n_evs_in_bldg = cars_per_building(bldg_id);
%     title(sprintf('Building %d | %d units | %d EVs', bldg_id, n_units, n_evs_in_bldg), ...
%          'FontWeight', 'bold');
%     
%     set(ax);
% end
% 
% 
% sgtitle('Attached Buildings: Load Profile & EV Discharge Validation', 'FontSize', 14, 'FontWeight', 'bold');
% % 
% %%
end