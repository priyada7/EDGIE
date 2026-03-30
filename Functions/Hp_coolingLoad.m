function [Hp_sizeCooling]=Hp_coolingLoad(thetaFull,solarFull,fullP,n1,Rlinearfit,designTempCool,floorArea,check)
% this function performs the simulation of heat pump with resistance heating
% equipments 
%



t1=1;
K=length(thetaFull);
R = Rlinearfit;
                       % conditioned floor area, m^2


thetaLow = f2c(82);                                      % first temperature point for heat pump COP curve, C
thetaHigh = f2c(95);                                    % second temperature point, C


etaLow = trirnd(4.00,4.40,1,n1);                            % first COP point
etaHigh = trirnd(3.2,3.30,1,n1);                 % second COP point
eta2 = -(etaLow + (etaHigh-etaLow).*(f2c(designTempCool)-thetaLow)./(thetaHigh-thetaLow)); % COP curve
%eta1(eta1 < 1)=1;


electricLoadThermalPower = fullP/n1; % thermal power from electrical loads other than heat pump, kW
bodyThermalPower = trirnd(0.1,0.5,K,t1); % thermal power from body heat, kW
solarThermalPower = 0.03*trirnd(0.9,1.1,K,t1).*floorArea.*solarFull; % thermal power from sunlight, kW 

    qe1 = (electricLoadThermalPower + ...
    bodyThermalPower + ...
    solarThermalPower);

if check == 1
   qe1=qe1./4;
end 

differences = abs(f2c(designTempCool) - thetaFull);
% Find the index of the minimum absolute difference
[~, ind] = min(differences);

Tset = f2c(trirnd(61,76,1,1));
thermalPower1= ( (Tset - f2c(designTempCool))/R - qe1(ind) );
Hp_sizeCooling = ((trirnd(1.2,1.3,1,1)*((abs(thermalPower1)./(0.80)))));

end