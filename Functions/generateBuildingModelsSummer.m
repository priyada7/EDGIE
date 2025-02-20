function [a1,C,R,eta1,Tset,qe,pMax1,pMaxAux]=generateBuildingModelsSummer(n1,dt,K,theta,t,P,solar,Rlinearfit,hpSize,floorArea)

% this function is used to generate input parameters for building model
% simulation
%
% Input:
%  n1, number of homes
%  dt, time sep, h
%  K, number of time steps
%  theta, Kx1 vector of outdoor temperature, C
%  t, (K+1)x1 vector of time span, h
%  P, Kx1 vector of base electricalload, kW
%  solar, Kx1 vector of solar irraditaion, kW/m^2
%  Rlinearfit, 1xn1 vector of thermal resistances, C/kW
%  floorArea, floor area, m^2
% hp size
% Output:
%  a1, 1xn1 matrix of discrete-time dynamics parameters
%  C, 1xn1 matrix of thermal capacitance, kWh/C
%  R, 1xn1 vector of thermal resistances, C/kW
%  eta1, Kxn1 matrix of heat pump COP
%  Tset, (K+1)xn1 matrix of heating temperature setpoint, C
%  qe, (K+1)xn1 matrix of exogenous thermal power, kW
%  pMax1, kxn1 matrix of max electical capcity of heat pump, kW
%  pMaxAux, kxn1 matrix of max resistance power, kW


%Building data
                   
ceilingHeight = repmat(3,1,n1);                   % ceiling height, m
rhoc = 3.42e-4;                                         % volumetric thermal capacitance of air, kWh/C/m^3
C = trirnd(16,18,1,n1).*(rhoc*floorArea'.*ceilingHeight) ;% thermal capacitance, kWh/C



R = Rlinearfit';

a1 = exp(-dt./(R.*C));                                  % discrete-time dynamics parameter


thetaLow = f2c(82);                                      % first temperature point for heat pump COP curve, C
thetaHigh = f2c(95); 


etaLow = trirnd(4.00,4.40,1,n1);                            % first COP point
etaHigh = trirnd(3.2,3.30,1,n1);                 % second COP point


eta1 = (etaLow + (etaHigh-etaLow).*(theta-thetaLow)./(thetaHigh-thetaLow)) ;% COP curve






%%%%%%%%%%%%%%%%%%%%%

Tset = repmat(f2c(trirnd(65,74,1,n1)),K+1,1);                             % heating temperature setpoint, C, tuned to cold-climate occupied values in RECS 2015 Table HC6.6

nSetback = 0;                                                             %round(n1/2); % number of units with night setbacks
for i=1:nSetback
    tDown = randi([22,23]);                             % setback time, h
    tUp = randi([7 10]);                                % recovery time, h
    setback = 5*trirnd(3,6,1,1)/9;                      % setback depth, C
    Tset(:,i) = Tset(:,i) - setback;                    % start with setback over all hours
    Tset(mod(t,24)>=tUp & mod(t,24)<=tDown,i) = Tset(mod(t,24)>=tUp & mod(t,24)<=tDown,i) + setback;
end


electricLoadThermalPower = P/n1; % thermal power from electrical loads other than heat pump, kW
bodyThermalPower = trirnd(0.1,0.5,K,n1); % thermal power from body heat, kW
solarThermalPower = 0.03*trirnd(0.9,1.1,K,n1).*floorArea'.*solar; % thermal power from sunlight, kW 
qe = electricLoadThermalPower + bodyThermalPower + solarThermalPower; % exogenous thermal power, kW (8-12 W/m^2)


pMax1 = hpSize.*ones(K,n1);                               % max electical capcity of heat pump
pMaxAux = 19.5.*ones(K,n1);                            % max electical capcity of resistance backup element

% fprintf('Mean exogenous thermal power intensity: %.3g W/m^2.\n',...
%    mean(mean(qe)./floorArea)*1000)
end