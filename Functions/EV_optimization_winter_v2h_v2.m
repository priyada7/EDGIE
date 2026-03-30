function [e,p2v,v2hEligible,HousingUnitsload,cvx_status]=...
    EV_optimization_winter_v2h_v2(K,p1baseHPWinter,p1basehpwhWinter,a2,e0,eMax,eMin,tau,ec,dt,pdMax,pcMax,atWork,atHome,...
    n2,Pwinter,pWorkBase,etac,etad,t,tw,th,onRoad,n1,winterPeak,TodaysHeadroom,FutureHeadroom,winterPeakComStock,U_units2bldg,A_attached_bldg,A_detached_home,PdetachedWinter,PattachedWinter,n1d,n1a,cars_attached,cars_detached,ev_start_attached,EnergyCommercial,EnergyResidential,normPrice,...
    EV_to_home,size_detached,size_attached,n_buildings,building_to_home,cars_per_building,frac,A_hv,C_bh,D_bv,scenario)
%--------------------------------------------------------------------------
% INPUT ARGUMENTS
%--------------------------------------------------------------------------
% K : scalar (integer)
%     Number of time steps (e.g., hours in the optimization horizon).
%
% p1baseHPSummer : K × n1 matrix
%     Baseline heat‑pump electric demand for each home.
%
% p1basehpwhSummer : K × n1 matrix
%     Baseline heat‑pump water‑heater electric demand for each home.
%
% a2 : n2 × 1 vector
%     EV battery self‑discharge/retention coefficient.
%
% e0 : 1 × n2 vector
%     Initial EV state of charge (SOC).
%
% eMax, eMin : 1 × n2 vectors
%     Maximum and minimum allowable EV SOC.
%
% tau : n2 × 1 vector
%     Charging efficiency scaling parameter in battery dynamics.
%
% ec : n2 × 1 vector
%     Energy consumed while driving (kWh).
%
% dt : scalar
%     Length of each time step (hours).
%
% pdMax, pcMax : n2 × 1 vectors
%     Maximum EV discharging and charging power limits (kW).
%
% atWork, atHome : K × n2 logical matrices
%     EV location indicators (1 if EV at work/home at time k).
%
% onRoad : K × n2 logical matrix
%     Indicator that EV is driving (no charge or discharge allowed).
%
% n2 : scalar
%     Total number of EVs.
%
% PSummer : K × 1 vector
%     Baseline residential non‑HVAC summer load (MELs).
%
% pWorkBase : K × 1 vector
%     Baseline commercial building load (ComStock).
%
% etac, etad : n2 × 1 vectors
%     EV charging and discharging efficiencies.
%
% t, tw, th : vectors
%     Time indexing information (not explicitly constrained here).
%
% n1 : scalar
%     Total number of homes.
%
% n1d, n1a : scalars
%     Number of detached and attached homes, respectively.
%
% PdetachedSummer : K × n1d matrix
%     Baseline summer load for detached homes.
%
% PattachedSummer : K × n1a matrix
%     Baseline summer load for attached units.
%
% U_units2bldg : n_buildings × n1a sparse matrix
%     Unit‑to‑building aggregation matrix for attached housing.
%
% A_attached_bldg : n_buildings × n2 sparse matrix
%     EV‑to‑attached‑building incidence matrix.
%
% A_detached_home : n1d × n2 sparse matrix
%     EV‑to‑detached‑home incidence matrix.
%
% cars_attached : n1a × 1 vector
%     EV count per attached unit.
%
% cars_detached : n1d × 1 vector
%     EV count per detached home.
%
% ev_start_attached : scalar
%     Index of first EV associated with attached buildings.
%
% SummerPeak : scalar
%     Reference summer peak load (kW).
%
% TodaysHeadroom, FutureHeadroom : scalars
%     Scaling factors for peak constraint enforcement.
%
% SummerPeakComStock : scalar
%     Commercial (workplace) summer peak load reference.
%
% EnergyResidential, EnergyCommercial : scalars
%     Energy cost weights for residential and commercial consumption.
%
% normPrice : scalar
%     Price coefficient for power ramping (l1‑norm penalty).
%
% v2hEligible : 1 × n2 logical vector
%     Indicates which EVs are allowed to discharge (V2H‑capable).
%
% A_hv : n1 × n2 sparse matrix
%     Home‑to‑EV incidence matrix.
%
% C_bh : (n1d+n_buildings) × n1 sparse matrix
%     Mapping from homes to detached homes and attached buildings.
%
% D_bv : (n1d+n_buildings) × n2 sparse matrix
%     Aggregated building/home‑to‑EV mapping.
%
% scenario : scalar {1,2,3}
%     V2H aggregation scenario:
%       1 → Building‑level V2H for attached homes
%       2 → Unit‑level V2H (no sharing)
%       3 → Detached‑only V2H (attached charging only)
%
%--------------------------------------------------------------------------
% OUTPUT ARGUMENTS
%--------------------------------------------------------------------------
% e : (K+1) × n2 matrix
%     Optimal EV state of charge trajectory.
%
% p2v : K × n2 matrix
%     Resulting EV power seen by buildings/homes:
%       • Positive → charging
%       • Negative → discharging (V2H)
%
% HousingUnitsload : K × n1 matrix
%     Baseline (non‑EV) electricity demand for all homes.
%
% cvx_status : string
%     Solver status returned by CVX (e.g., 'Solved', 'Infeasible').

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



% --- split EV indices ---
ev_det_idx = 1:(ev_start_attached-1);        % detached EVs
ev_att_idx = ev_start_attached:n2;           % attached EVs

if scenario == 3
pChemMin(:, ev_att_idx) = 0;
end
%%

% --- Disable V2H for (1 - frac) of EVs ---
    
    mEnable  = round(frac * n2);         % EVs that CAN do V2H
    idxPerm  = randperm(n2);
    idxEnable = idxPerm(1:mEnable);
    v2hEligible = false(1, n2);
    v2hEligible(idxEnable) = true;

    % Set non-eligible EVs to no-discharge
    pChemMin(:, ~v2hEligible) = 0;





%%
% When EV is at work: only charging allowed (no V2H at work)
pChemMax = pChemMax + atWork .* repmat(valCharge, K, 1);

% When EV is on road: no charging/discharging (set limits to zero)
pChemMin(onRoad) = 0;
pChemMax(onRoad) = 0;

% Calculate energy consumption while driving
for k=1:K
   w(k,onRoad(k,:)) = ec(onRoad(k,:))/dt;
   %pChemMin(k,onRoad(k,:)) = -etad(onRoad(k,:)).*ec(onRoad(k,:))/dt;
end

% --------------------------
% --------------------------
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
fprintf('Setting up CVX optimization for winter...\n');

%%
cvx_begin quiet
cvx_precision low

    variables e(K+1,n2) pChem(K,n2)

expressions homeLoadWinter workLoadWinter diffCalc;  % just declare holder first

% Total loads
homeLoadWinter = Pwinter...% MEL MFRED
            + sum(p1baseHPWinter,2)... % Heat pump
            + sum(p1basehpwhWinter,2) ...% Heat-pump water heater
            + sum(atHome.*max(etad.*pChem, pChem./etac),2); % EV

workLoadWinter = pWorkBase + sum(atWork .* (pChem ./ etac), 2); % ComStock + EVcharging
diffCalc=pChem(2:K,:) - pChem(1:K-1,:);

minimize( ...
    DeviationPenalty_EV*sum(max(0, e0 - e(K+1,:))) ...   % EV deviation penalty
    + (1/K   )*EnergyResidential*sum(homeLoadWinter) ...                          % Energy usage cost
    + (1/K)*EnergyCommercial*sum(workLoadWinter) ...                          % Energy usage cost
    + 1000*pos(FutureHeadroom*max(homeLoadWinter) - TodaysHeadroom*winterPeak) ...            % Winter peak violation penalty
    + 1000*pos(FutureHeadroom*max(workLoadWinter) - TodaysHeadroom*winterPeakComStock) ... 
    + (1/(K-1))*normPrice*norm(diffCalc, 1) ... 
)  


  


    subject to
    % --------------------------
    % Constraints
    % --------------------------
    e(1,:) == e0;
    eMin <= e <= eMax;
    e(2:K+1,:) == a2.*e(1:K,:) + (1-a2).*tau.*(pChem - w);


    pChemMin <= pChem <= pChemMax;

    % DETACHED: Per-home demand (K × n_detached)
    homeDemand_detached = PdetachedWinter+ ...
                         p1baseHPWinter(:,1:n_detached) + ...
                         p1basehpwhWinter(:,1:n_detached);
    
    % ATTACHED: Per-building demand (K × n_buildings)
    % Aggregate all unit demands within each building
    homeDemand_attached_units =PattachedWinter+ ...
                                p1baseHPWinter(:,n_detached+1:n1) + ...
                               p1basehpwhWinter(:,n_detached+1:n1);  % K × n1a
    
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
    elseif scenario == 4
    % V2G: no home load constraint, all EVs can discharge
    % already handled by pChemMin bounds — just skip constraint
end

cvx_end
e=full(e);
pChem=full(pChem);
p2v = max(etad .* pChem, pChem ./ etac) .* atHome ...   % V2H: at home
    + (pChem ./ etac) .* atWork;                          % charging: at work
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
%     home_load = PdetachedWinter(:, home_id) + p1baseHPWinter(:, home_id) + p1basehpwhWinter(:, home_id);
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
% 
% %%
end