function [Hp_sizeCooling]=Hp_coolingLoad(thetaFull,solarFull,fullP,n1,Rlinearfit,designTempCool,floorArea)
% this function is used to size the heat pump for cooling
% equipments 
%
% Input:
%  thetaFull, Kx1 full year outdoor temp, C
%  solarFull, Kx1 full year solar radiation
%  fullP, Kx1 misc load
%  n1, number of homes
%  Rlinearfit, 1xn1 vector of thermal resistances, C/kW
%  designTempCool, cooling design temp, C
%  floorArea, floor area, m^2
%
% Output:
%  Hp_sizeCooling, electric capacity of heat pump, kW



t1=1;
K=length(thetaFull);
R = Rlinearfit;
                       % conditioned floor area, m^2


thetaLow = f2c(82);                                      % first temperature point for heat pump COP curve, C
thetaHigh = f2c(95);                                    % second temperature point, C
etaLow = repmat(4.24,1,t1);                            % first COP point
etaHigh = repmat(3.5,1,t1); 
eta2 = -(etaLow + (etaHigh-etaLow).*(f2c(designTempCool)-thetaLow)./(thetaHigh-thetaLow)); % COP curve


electricLoadThermalPower = fullP/n1; % thermal power from electrical loads other than heat pump, kW
bodyThermalPower = gauss(0.1,0.5,K,t1); % thermal power from body heat, kW
solarThermalPower = 0.03*gauss(0.9,1.1,K,t1).*floorArea.*solarFull; % thermal power from sunlight, kW 
qe1 = electricLoadThermalPower + ...
    bodyThermalPower + ...
    solarThermalPower;
differences = abs(f2c(designTempCool) - thetaFull);
% Find the index of the minimum absolute difference
[~, ind] = min(differences);

Tset = f2c(trirnd(61,76,1,1));
thermalPower1= ( (Tset - f2c(designTempCool))/R - qe1(ind) );

Hp_sizeCooling = ceil(trirnd(1.2,1.4,1,1)*ceil((thermalPower1./(eta2*0.85))));

end