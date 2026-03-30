function [U_units2bldg,A_attached_bldg,A_detached_home,cars_attached,cars_detached,ev_start_attached,EV_to_home,size_detached,size_attached,n_buildings,building_to_home,cars_per_building,...
    A_hv,C_bh,D_bv] = EV_Allocation_v2(units_per_building,n1,n2,n1a,n1d,floorArea,ft2m2)

% % Complete EV allocation with realistic probability-based distribution
% and building-level sharing for attached homes


%--------------------------------------------------------------------------
% INPUT ARGUMENTS
%--------------------------------------------------------------------------
% units_per_building : scalar (positive integer) Average number of attached units per building.
%
% n1 : Total number of homes (detached + attached).
%
% n2 : Total number of EVs to be allocated across all homes.
%
% n1a : Number of attached homes (e.g., apartments, townhomes).
%
% n1d : Number of detached homes (single‑family).
%
% floorArea : 1 × n1 vector Floor area of each home 
%
% ft2m2 : Conversion factor from ft² to m² (e.g., 10.7639).
%
%--------------------------------------------------------------------------
% OUTPUT ARGUMENTS
%--------------------------------------------------------------------------
% cars_detached : n1d × 1 vector
%     Number of EVs assigned to each detached home.
%
% cars_attached : n1a × 1 vector
%     Number of EVs assigned to each attached unit.
%
% size_detached : n1d × 1 vector
%     Floor area of detached homes (m²).
%
% size_attached : n1a × 1 vector
%     Floor area of attached units (m²).
%
% n_buildings : scalar (integer)
%     Total number of attached buildings created from attached units.
%
% building_to_home : n_buildings × 1 cell array
%     Each cell contains global home indices belonging to that building.
%
% cars_per_building : n_buildings × 1 vector
%     Total number of EVs assigned to each attached building.
%
% EV_to_home : 1 × n2 vector
%     Mapping of each EV to:
%       • detached home index (1 … n1d), or
%       • attached building index (1 … n_buildings)
%
% ev_start_attached : scalar
%     Index of first EV assigned to attached buildings in EV_to_home.
%
%--------------------------------------------------------------------------
% INCIDENCE / AGGREGATION MATRICES (SPARSE)
%--------------------------------------------------------------------------
% A_detached_home : n1d × n2 sparse matrix
%     Detached‑home‑to‑EV mapping (one EV per column).
%
% A_attached_bldg : n_buildings × n2 sparse matrix
%     Building‑level mapping of EVs in attached housing.
%
% U_units2bldg : n_buildings × n1a sparse matrix
%     Maps attached units to their corresponding buildings.
%
% A_hv : n1 × n2 sparse matrix
%     Unified home‑to‑EV incidence matrix:
%       rows → homes
%       cols → EVs
%
% C_bh : (n1d + n_buildings) × n1 sparse matrix
%     Maps homes to aggregation nodes:
%       rows 1…n1d              → detached homes (identity)
%       rows n1d+1…end          → attached buildings
%
% D_bv : (n1d + n_buildings) × n2 sparse matrix
%     Final building‑/home‑to‑EV aggregation:
%         D_bv = C_bh * A_hv
%

    
% ------------------ HOME CLASSIFICATION ------------------
n_detached = n1d;
n_attached = n1a ;
n_buildings = ceil(n_attached / units_per_building);

% fprintf('Total homes: %d\n', n1);
% fprintf('Detached homes: %d (homes 1-%d)\n', n_detached, n_detached);
% fprintf('Attached units: %d (homes %d-%d)\n', n_attached, n_detached+1, n1);
% fprintf('Attached buildings: %d (avg %.1f units/building)\n', n_buildings, n_attached/n_buildings);
% fprintf('Total EVs available: %d\n', n2);


% ------------------ GENERATE HOME SIZES ------------------

size_detached = floorArea(1,1:n_detached)'./ft2m2;  % 1000–5500 sq ft
size_attached = floorArea(1,n_detached+1:end)'./ft2m2;  % 1500–3000 sq ft

% ------------------ ALLOCATE CARS TO DETACHED HOMES ------------------
%fprintf('\n--- Phase 1: Allocating to detached homes (probability-based) ---\n');

% Probability table from research data
prob_table = [...
    1200, [0.70 0.18 0.08 0.03];  % <1200 sqft
    1600, [0.30 0.17 0.51 0.02];  % 1200-1600
    2000, [0.10 0.12 0.72 0.06];  % 1600-2000
    2400, [0.05 0.06 0.80 0.09];  % 2000-2400
    3000, [0.02 0.02 0.78 0.18];  % 2400-3000
    5000, [0.01 0.03 0.59 0.37];  % 3000-5000
    6000, [0.01 0.02 0.30 0.67]]; % >5000

cars_detached = zeros(n_detached,1);
for i = 1:n_detached
    sqft = size_detached(i);
    % Find matching range
    idx = find(sqft <= prob_table(:,1),1,'first');
    if isempty(idx), idx = size(prob_table,1); end
    
    p = prob_table(idx,2:5);  % Probabilities for 0,1,2,3 cars
    outcomes = [0 1 2 3];
    cars_detached(i) = randsample(outcomes,1,true,p);
end

cars_assigned = sum(cars_detached);
if cars_assigned > n2
    %fprintf('NOTE: sampled %d cars for detached homes but only %d total EVs available. Trimming detached assignments...\n', cars_assigned, n2);
    excess = cars_assigned - n2;
    % repeatedly remove 1 car from a random detached home that has >0 cars until excess removed
    idx_nonzero = find(cars_detached>0);
    while excess > 0 && ~isempty(idx_nonzero)
        sel = idx_nonzero(randi(length(idx_nonzero)));
        cars_detached(sel) = cars_detached(sel) - 1;
        excess = excess - 1;
        idx_nonzero = find(cars_detached>0); % update
    end
    cars_assigned = sum(cars_detached);
    %fprintf('After trimming: cars_assigned (detached) = %d\n', cars_assigned);
end
 cars_remaining = n2 - cars_assigned;
if n_attached == 0 && cars_remaining > 0
   % fprintf('No attached homes. Allocating %d remaining cars to detached homes by lowest cars/sqft ratio...\n', cars_remaining);
    
    while cars_remaining > 0
        % Calculate cars per square foot ratio for all detached homes
        ratios = cars_detached ./ size_detached;
        
        % Find home with minimum ratio (lowest cars per sqft)
        [~, idx_min] = min(ratios);
        
        % Allocate one car to that home
        cars_detached(idx_min) = cars_detached(idx_min) + 1;
        cars_remaining = cars_remaining - 1;
    end
end
%     fprintf('All cars allocated to detached homes: %d total\n', sum(cars_detached));
% fprintf('Cars assigned to detached homes: %d\n', cars_assigned);
% fprintf('Average cars per detached home: %.2f\n', mean(cars_detached));




% ------------------ ALLOCATE CARS TO ATTACHED UNITS  ------------------
%fprintf('\n--- Phase 2: Allocating to attached units (floor area-based) ---\n');

cars_remaining = n2 - cars_assigned;
%fprintf('Cars remaining for attached units: %d\n', cars_remaining);

% Allocate EVs to each attached unit based on floor area (like detached, but simpler)
cars_attached = zeros(n_attached, 1);

for i = 1:n_attached
    sqft = size_attached(i);
    
    % Probability of EV ownership for attached units (simpler than detached)
    if sqft < 1000
        % Small attached units: mostly 0-1 EVs
        p = [0.60, 0.40, 0.00, 0.00];  % 60% no EV, 35% one EV, 5% two EVs
    else
        % Large attached units
        p = [0.15, 0.50, 0.35, 0.00];  % More likely to have 1-2 EVs
    end
    outcomes = [0 1 2 3];
    cars_attached(i) = randsample(outcomes, 1, true, p);
end


  if n_attached>0  
    % Distribute remainder
    deficit = cars_remaining - sum(cars_attached);
    for i = 1:deficit
        % Give to units that originally had EVs
        idx = find(cars_attached < 2);
        if ~isempty(idx)
            selected = idx(randi(length(idx)));
            cars_attached(selected) = cars_attached(selected) + 1;
        end
    end

     still_deficit = cars_remaining - sum(cars_attached);
     while still_deficit>0
         ratios = cars_attached ./ size_attached;
         [~, idx_min] = min(ratios);

         % Allocate one car to that home
         cars_attached(idx_min) = cars_attached(idx_min) + 1;

         % Update deficit
         still_deficit = still_deficit - 1;
     end
  end
total_attached_sampled = sum(cars_attached);
if total_attached_sampled > cars_remaining
    excess = total_attached_sampled - cars_remaining;
 %   fprintf('NOTE: sampled %d attached EVs but only %d cars remain. Trimming attached assignments...\n', total_attached_sampled, cars_remaining);
    idx_nonzero = find(cars_attached>0);
    while excess > 0 && ~isempty(idx_nonzero)
        sel = idx_nonzero(randi(length(idx_nonzero)));
        cars_attached(sel) = cars_attached(sel) - 1;
        excess = excess - 1;
        idx_nonzero = find(cars_attached>0);
    end
end



% ------------------ NOW GROUP UNITS INTO BUILDINGS ------------------
%fprintf('\n--- Grouping attached units into buildings ---\n');

% Create building assignment
building_assignment = zeros(n_attached, 1);
for i = 1:n_attached
    building_assignment(i) = ceil(i / units_per_building);
end
n_buildings = max(building_assignment);
n_buildings = max([0; building_assignment]);   % always returns scalar 0 if empty
% Create building_to_home mapping
building_to_home = cell(n_buildings, 1);
for b = 1:n_buildings
    units_in_building = find(building_assignment == b);
    % Convert to actual home indices (offset by n_detached)
    building_to_home{b} = n_detached + units_in_building';
end

% Calculate EVs per building (sum of all units in that building)
cars_per_building = zeros(n_buildings, 1);
for b = 1:n_buildings
    units_in_building = find(building_assignment == b);
    cars_per_building(b) = sum(cars_attached(units_in_building));
end




EV_to_home = zeros(1, n2);
EV_count = 0;

% Phase 1: Map EVs to detached homes (one-to-one)
for home = 1:n_detached
    num_cars = cars_detached(home);
    for car = 1:num_cars
        if EV_count < n2
            EV_count = EV_count + 1;
            EV_to_home(EV_count) = home;
        end
    end
end
%fprintf('EVs mapped to detached homes: %d\n', EV_count);

% Phase 2: Map EVs to buildings
% Note: We assign EVs based on individual unit ownership, but map to building
ev_start_attached = EV_count + 1;

for b = 1:n_buildings
    num_cars = cars_per_building(b);  % Total EVs in this building
    for car = 1:num_cars
        if EV_count < n2
            EV_count = EV_count + 1;
            EV_to_home(EV_count) = b;  % Store building ID
        end
    end
end


% Handle unmapped EVs
if EV_count < n2
   % fprintf('WARNING: %d EVs unmapped. Distributing to buildings...\n', n2 - EV_count);
    b = 1;
    for ev = (EV_count+1):n2
        EV_to_home(ev) = b;
        b = b + 1;
        if b > n_buildings
            b = 1;
        end
    end
    EV_count = n2;
end



% ---------- DETACHED: homes × EVs (one-to-one) ----------
row_d = [];
col_d = [];
val_d = [];

for ev = 1:ev_start_attached-1
    home = EV_to_home(ev);         % detached home id in 1..n_detached
    row_d(end+1) = home;
    col_d(end+1) = ev;
    val_d(end+1) = 1;
end

A_detached_home = sparse(row_d, col_d, val_d, n_detached, EV_count);%n1dxn2
%fprintf('A_detached_home: %d × %d, nnz=%d\n', size(A_detached_home,1), size(A_detached_home,2), nnz(A_detached_home));

% ---------- ATTACHED: buildings × EVs (EVs mapped to their building) ----------
if n_attached >0

    row_b = [];
    col_b = [];
    val_b = [];

    for ev = ev_start_attached:EV_count
        bldg = EV_to_home(ev); % building id in 1..n_buildings
        if 1 <= bldg && bldg <= n_buildings
            row_b(end+1) = bldg;
            col_b(end+1) = ev;
            val_b(end+1) = 1;
        end
    end

    A_attached_bldg = sparse(row_b, col_b, val_b, n_buildings, EV_count);%bldgxn2
 %   fprintf('A_attached_bldg: %d × %d, nnz=%d\n', size(A_attached_bldg,1), size(A_attached_bldg,2), nnz(A_attached_bldg));

    % ---------- ATTACHED aggregator: buildings × attached-units ----------
    % Build units→buildings matrix to aggregate unit-level demand into building totals.
    n_attached_units = n1 - n_detached;

    row_u = [];
    col_u = [];
    val_u = [];

    for b = 1:n_buildings
        homes = building_to_home{b};                   % global home indices
        homes = homes(homes >= n_detached+1 & homes <= n1);  % ensure attached only
        if ~isempty(homes)
            u_idx = homes - n_detached;               % convert to 1..n_attached_units
            row_u = [row_u, repmat(b, 1, numel(u_idx))];
            col_u = [col_u, u_idx];
            val_u = [val_u, ones(1, numel(u_idx))];
        end
    end

    U_units2bldg = sparse(row_u, col_u, val_u, n_buildings, n_attached_units);%bldgxn1a
  %  fprintf('U_attached_bldg: %d × %d, nnz=%d\n', size(U_units2bldg,1), size(U_units2bldg,2), nnz(U_units2bldg));
else
    A_attached_bldg = zeros;
    U_units2bldg = zeros;

end

% ------------------ BUILD UNIFIED HOME×EV MATRIX ------------------
% A_hv: n1 × n2 (each row = home, each col = EV)

% --- DETACHED block (rows 1..n_detached) ---
% A_detached_home already has correct rows 1..n_detached vs all n2 EVs
A_detached_block = A_detached_home;  % n_detached × n2

% --- ATTACHED block (rows n_detached+1 .. n1) ---
% We need to map each EV to its specific attached UNIT, not just building.
% Strategy: within each building, assign EVs round-robin (or by unit ownership)
% using cars_attached to know how many EVs each unit owns.

row_a = [];
col_a = [];
val_a = [];

for b = 1:n_buildings
    % Find which EVs belong to this building
    ev_in_bldg = find(EV_to_home(ev_start_attached:end) == b) + ev_start_attached - 1;
    
    % Find which units are in this building (local indices into cars_attached)
    units_in_b = find(building_assignment == b);  % local 1..n_attached indices
    
    % Assign EVs to units based on cars_attached ownership
    ev_ptr = 1;
    for j = 1:length(units_in_b)
        unit_local = units_in_b(j);
        unit_global = n_detached + unit_local;  % global home index
        n_ev_unit = cars_attached(unit_local);
        
        for k = 1:n_ev_unit
            if ev_ptr <= length(ev_in_bldg)
                row_a(end+1) = unit_global;
                col_a(end+1) = ev_in_bldg(ev_ptr);
                val_a(end+1) = 1;
                ev_ptr = ev_ptr + 1;
            end
        end
    end
end

A_attached_unit = sparse(row_a, col_a, val_a, n1, n2);  % n1 × n2 (only attached rows filled)

% --- COMBINE: stack detached block on top ---
% Detached block already spans all n2 columns
A_detached_full = sparse(double(row_d), double(col_d), double(val_d), n1, n2);  % pad to n1 rows

A_hv = A_detached_full + A_attached_unit;  % n1 × n2

% ------------------ BUILD C MATRIX (buildings × homes) ------------------
% C: (n_detached + n_buildings) × n1
% Rows 1..n_detached          → detached home i maps to home i (identity)
% Rows n_detached+1..n_detached+n_buildings → attached building maps to its units

row_c = [];
col_c = [];
val_c = [];

% --- DETACHED: identity block ---
for i = 1:n_detached
    row_c(end+1) = i;
    col_c(end+1) = i;
    val_c(end+1) = 1;
end

% --- ATTACHED: each building points to its unit home indices ---
for b = 1:n_buildings
    homes = building_to_home{b};   % global home indices (n_detached+1 .. n1)
    for h = homes
        row_c(end+1) = n_detached + b;
        col_c(end+1) = h;
        val_c(end+1) = 1;
    end
end

C_bh = sparse(row_c, col_c, val_c, n_detached + n_buildings, n1);


D_bv = C_bh*A_hv;

end