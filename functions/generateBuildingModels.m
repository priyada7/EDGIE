function [a1,C,R,eta1,Tset,qe,pMax1,pMaxAux]=generateBuildingModels(n1,dt,K,theta,t,P,solar)

%Building data
floorArea = trirnd(200,250,1,n1);                       % conditioned floor area, m^2
ceilingHeight = trirnd(3,3.6,1,n1);                     % ceiling height, m
rhoc = 3.42e-4;                                         % volumetric thermal capacitance of air, kWh/C/m^3
C = trirnd(16,18,1,n1).*(rhoc*floorArea.*ceilingHeight); % thermal capacitance, kWh/C

thetaBase = f2c(50);                                    % base temperature for heading degree-hours, C
RValue = gauss(550,550,1,n1);                                           % m^2*C/kW, value obtained from DC house validation
R = RValue./(4.*ceilingHeight.*sqrt(2.*floorArea));     % thermal resistance, C/kW

a1 = exp(-dt./(R.*C));                                  % discrete-time dynamics parameter
thetaLow = f2c(5);                                      % first temperature point for heat pump COP curve, C
thetaHigh = f2c(45);                                    % second temperature point, C
etaLow = trirnd(1.5,2,1,n1);                            % first COP point
etaHigh = etaLow + 1.5 + trirnd(0,0.5,1,n1);            % second COP point
eta1 = etaLow + (etaHigh-etaLow).*(theta-thetaLow)./(thetaHigh-thetaLow); % COP curve
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
bodyThermalPower = gauss(0.1,0.5,K,n1); % thermal power from body heat, kW
solarThermalPower = 0.03*gauss(0.9,1.1,K,n1).*floorArea.*solar; % thermal power from sunlight, kW 
qe = electricLoadThermalPower + bodyThermalPower + solarThermalPower; % exogenous thermal power, kW (8-12 W/m^2)


pMax1 = 4.5*ones(K,n1);                               % max electical capcity of heat pump
pMaxAux = 19.5*ones(K,n1);                            % max electical capcity of resistance backup element

fprintf('Mean exogenous thermal power intensity: %.3g W/m^2.\n',...
    mean(mean(qe)./floorArea)*1000)
end