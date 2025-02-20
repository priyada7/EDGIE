function [Hp_sizeHeating]=Hp_heatingLoad(thetaFull,solarFull,fullP,n1,Rlinearfit,designTempHeat,floorArea,zone)
% this function is used to size the heat pump for heating
% equipments 
%
% Input:
%  thetaFull, Kx1 full year outdoor temp, C
%  solarFull, Kx1 full year solar radiation
%  fullP, Kx1 misc load
%  n1, number of homes
%  Rlinearfit, 1xn1 vector of thermal resistances, C/kW
%  designTempHear, heating design temp, C
%  floorArea, floor area, m^2
%  zone, ASHRAE zone location
% Output:
%  Hp_sizeCooling, electric capacity of heat pump, kW


t1=1;

K=length(thetaFull);
R = Rlinearfit;
                      % conditioned floor area, m^2

thetaLow = f2c(5);                                      % first temperature point for heat pump COP curve, C
thetaHigh = f2c(45);                                    % second temperature point, C

if zone > 4
        etaLow = trirnd(2.2,2.4,1,t1);                            % first COP point
        etaHigh = etaLow + 1.5 + trirnd(0,0.5,1,t1);            % second COP point
        eta2 = etaLow + (etaHigh-etaLow).*(f2c(designTempHeat)-thetaLow)./(thetaHigh-thetaLow); % COP curve
        eta2(eta2 < 1)=1;
    else
        etaLow = trirnd(1.5,2,1,t1);                            % first COP point
        etaHigh = etaLow + 1.5 + trirnd(0,0.5,1,t1);            % second COP point
         eta2 = etaLow + (etaHigh-etaLow).*(f2c(designTempHeat)-thetaLow)./(thetaHigh-thetaLow); % COP curve
         eta2(eta2 < 1)=1;
end

electricLoadThermalPower = fullP/n1; % thermal power from electrical loads other than heat pump, kW
bodyThermalPower = gauss(0.1,0.5,K,t1); % thermal power from body heat, kW
solarThermalPower = 0.03*gauss(0.9,1.1,K,t1).*floorArea.*solarFull; % thermal power from sunlight, kW 
qe1 = electricLoadThermalPower + ...
    bodyThermalPower + ...
    solarThermalPower; % exogenous thermal power, kW (8-12 W/m^2)

[minTemp,idx]= min(thetaFull);

%Tset = cell2mat(Tsetcombined);
Tset = f2c(trirnd(61,76,1,1));
DesignTemp = f2c(designTempHeat);
 
thermalPower1= ( (Tset - DesignTemp)/R - qe1(idx) );
   

Hp_sizeHeating = ceil(trirnd(1.2,1.4,1,1)*(thermalPower1./eta2));

end