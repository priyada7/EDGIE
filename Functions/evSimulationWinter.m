function [phBase,pwBase,eBase,pdBase]=evSimulationWinter(K,n2,e0,p0,eMin,eMax,pcMax,a2,tau,etac,etad,dt,ec,onRoad,atHome,atWork)

% this function is used to run the simulation of electric vehicles 
%
% Input:
%  K, number of time steps
%  n2, number of ev
%  e0, initial energy of ev battery, kWh
%  eMin, n2xK matrix of min-energy capacity, kWh 
%  eMax, n2xK matrix of max-energy capacity, kWh 
%  pcMax, n2xK matrix of charge capacity, kW
%  a2, 1xn2 vector of discrete-time dynamics parameters
%  tau, KxL matrix of dissipation rate, 1/h
%  etac, 1xn2 vector of charge efficiency
%  etad, 1xn2 vector of discharge efficiency
%  dt, time step, h
%  ec, 1xn2 vector of commute energy, kWh
%  onRoad, indicator that vehicle's on the road
%  atHome,indicator that vehicle's at home
%  atWork, indicator that vehicle's at work
%
%
% Output:
%  phBase, Kxn2 matrix of power consumption from ev at home, kW
%  eBase, (K+1)xn2 matrix of stored energy in battery, kWh
%  pwBase, Kxn2 matrix of power consumption from ev at work, kW
%  pdBase, Kxn2 matrix of ev discharge power, kW

eBase = zeros(K+1,n2); % stored energy, kWh
eBase(1,:) = e0;
pcBase = zeros(K,n2); % home or work charge power, kW
pdBase = zeros(K,n2); % discharge power, kW
for k=1:K
    % initialize charging power with previous charging power
    if k==1
        pck = p0;
    else
        pck = pcBase(k-1,:);
    end
    
    % start charging if energy dropped below minimum
    pck(eBase(k,:)<eMin) = pcMax(eBase(k,:)<eMin);
    
    % stop charging if driving
    pck(~atHome(k,:) & ~atWork(k,:)) = 0;
    
    % saturate charging power
    pcBase(k,:) = max(0,min(pck,(eMax-a2.*eBase(k,:))./(1-a2)./tau./etac));
    
    % discharge power if on the road
    pdBase(k,onRoad(k,:)) = etad(onRoad(k,:)).*ec(onRoad(k,:))/dt;
    
    % update stored energy
    eBase(k+1,:) = a2.*eBase(k,:) + (1-a2).*tau.*(etac.*pcBase(k,:) - pdBase(k,:)./etad);
end
phBase = atHome.*pcBase; % baseline home charge power, kW
pwBase = atWork.*pcBase; % baseline work charge power, kW



end