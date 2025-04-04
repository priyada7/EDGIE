function [Hp_sizeHeating]=Hp_heatingLoad(thetaFull,solarFull,fullP,n1,Rlinearfit,designTempHeat,floorArea,zone,check)
% this function performs the simulation of heat pump with resistance heating
% equipments 
%
% Input:
%  K, number of time steps
%  n1, number of homes
%  pMax1, kxn1 matrix of max electical capcity of heat pump, kW
%  eta1, Kxn1 matrix of heat pump COP
%  pMaxAux, kxn1 matrix of max resistance power, kW
%  Tset, (K+1)xn1 matrix of heating temperature setpoint, C
%  qe, (K+1)xn1 matrix of exogenous thermal power, kW
%  a1, 1xn1 matrix of discrete-time dynamics parameters
%  R, 1xn1 vector of thermal resistances, C/kW
%  theta, Kx1 vector of outdoor temperature, C
%
%
% Output:
%  p1base, Kxn1 matrix for heat pump electric power, kW
%  HPload, n1x1 vector for heat pump electrical power, kW
%  Tbase2, (K+1)xn1 matrix for indoor temperature, C

% this function performs the simulation of heat pump with resistance heating
% equipments 
%
% Input:
%  K, number of time steps
%  n1, number of homes
%  pMax1, kxn1 matrix of max electical capcity of heat pump, kW
%  eta1, Kxn1 matrix of heat pump COP
%  pMaxAux, kxn1 matrix of max resistance power, kW
%  Tset, (K+1)xn1 matrix of heating temperature setpoint, C
%  qe, (K+1)xn1 matrix of exogenous thermal power, kW
%  a1, 1xn1 matrix of discrete-time dynamics parameters
%  R, 1xn1 vector of thermal resistances, C/kW
%  theta, Kx1 vector of outdoor temperature, C
%
%
% Output:
%  p1base, Kxn1 matrix for heat pump electric power, kW
%  HPload, n1x1 vector for heat pump electrical power, kW
%  Tbase2, (K+1)xn1 matrix for indoor temperature, C

% March - Sep end
t1=1;
K=length(thetaFull);
R = Rlinearfit;
              
electricLoadThermalPower = fullP/n1; % thermal power from electrical loads other than heat pump, kW
bodyThermalPower = trirnd(0.1,0.5,K,t1); % thermal power from body heat, kW
solarThermalPower = 0.03*trirnd(0.9,1.1,K,t1).*floorArea.*solarFull; % thermal power from sunlight, kW 
qe1 = electricLoadThermalPower + ...
    bodyThermalPower + ...
    solarThermalPower; % exogenous thermal power, kW (8-12 W/m^2)
if check == 1
   qe1=qe1./4;
end 
[minTemp,idx]= min(thetaFull);

%Tset = cell2mat(Tsetcombined);
Tset = f2c(trirnd(61,76,1,1));
DesignTemp = f2c(designTempHeat);
 
thermalPower1= ( (Tset - DesignTemp)/R - qe1(idx) );
Hp_sizeHeating = (trirnd(1.2,1.3,1,1)*(thermalPower1));



end