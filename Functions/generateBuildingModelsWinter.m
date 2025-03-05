function [eta1,qe]=generateBuildingModelsWinter(zone,n1,K,theta,P,solar,floorArea)

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
%
% Output:
%  a1, 1xn1 matrix of discrete-time dynamics parameters
%  C, 1xn1 matrix of thermal capacitance, kWh/C
%  R, 1xn1 vector of thermal resistances, C/kW
%  eta1, Kxn1 matrix of heat pump COP
%  Tset, (K+1)xn1 matrix of heating temperature setpoint, C
%  qe, (K+1)xn1 matrix of exogenous thermal power, kW
%  pMax1, kxn1 matrix of max electical capcity of heat pump, kW
%  pMaxAux, kxn1 matrix of max resistance power, kW

thetaLow = f2c(5);                                      % first temperature point for heat pump COP curve, C
thetaHigh = f2c(45); 

% etaLow = trirnd(3,4,1,n1);                            % first COP point
% etaHigh = etaLow + 1.5 + trirnd(0,0.5,1,n1);            % second COP point
%  eta1 = etaLow + (etaHigh-etaLow).*(theta-thetaLow)./(thetaHigh-thetaLow); % COP curve
%  eta1(eta1 < 1)=1;

    if zone > 4
       
        etaLow = trirnd(2.2,2.4,1,n1);                            % first COP point 2.2 2.4 trirnd(1.5,1.8,1
        etaHigh =etaLow + 1.5 + trirnd(0,0.5,1,n1);            % second COP point
        eta1 = etaLow + (etaHigh-etaLow).*(theta-thetaLow)./(thetaHigh-thetaLow); % COP curve
        eta1(eta1 < 1)=1;
    else
        etaLow =  trirnd(1.5,2.0,1,n1);                             % first COP point  trirnd(1.5,2,1,n1); trirnd(1.1,1.5,   1,homeWithHP);    
        etaHigh = etaLow + 1.5 + trirnd(0,0.5,1,n1);            % second COP point
        eta1 = etaLow + (etaHigh-etaLow).*(theta-thetaLow)./(thetaHigh-thetaLow); % COP curve
        eta1(eta1 < 1)=1;
    end


%%%%%%%%%%%%%%%%%%%%%

electricLoadThermalPower = P/n1; % thermal power from electrical loads other than heat pump, kW
bodyThermalPower = trirnd(0.1,0.5,K,n1); % thermal power from body heat, kW
solarThermalPower = 0.03*trirnd(0.9,1.1,K,n1).*floorArea'.*solar; % thermal power from sunlight, kW 
qe = electricLoadThermalPower + bodyThermalPower + solarThermalPower; % exogenous thermal power, kW (8-12 W/m^2)


end