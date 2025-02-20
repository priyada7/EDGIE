function [electricPower,summerPeak,waterHeatingLoad]=summerPeakcalc(thetaFull,solarFull,dt,fullP,n1,weatherTime,R,pMax1,HPWHloadFullyear,homeWithElectricWH,HPWHeachHouseLoad,floorArea)
% this function performs the simulation of heat pump with resistance
% heating,n1
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
t1=n1;
%R=R(1,1);

K=length(thetaFull);
Tbase2 = zeros(K+1,t1);      % indoor temperature, C
Tset = repmat(f2c(trirnd(65,74,1,t1)),K+1,1); 
Tbase2(1,:) = Tset(1,:);

                    % conditioned floor area, m^2
ceilingHeight = repmat(3,1,t1);
rhoc = 3.42e-4;                                         % volumetric thermal capacitance of air, kWh/C/m^3
C = trirnd(16,18,1,t1).*(rhoc*floorArea'.*ceilingHeight); % thermal capacitance, kWh/C
a1 = exp(-dt./(R.*C));
   
thetaLow = f2c(82);                                      % first temperature point for heat pump COP curve, C
thetaHigh = f2c(95);                                    % second temperature point, C

% etaLow = repmat(4.3,1,t1);                            % first COP point
% etaHigh = repmat(2.2,1,t1); 
% eta1 = (etaLow + (etaHigh-etaLow).*(thetaFull-thetaLow)./(thetaHigh-thetaLow)) ;% COP curve
      
etaLow = trirnd(3.5,4.00,1,n1);                            % first COP point
etaHigh = trirnd(2.5,2.8,1,n1);                 % second COP point


eta1 = (etaLow + (etaHigh-etaLow).*(thetaFull-thetaLow)./(thetaHigh-thetaLow)); % COP curve


electricLoadThermalPower = fullP/n1; % thermal power from electrical loads other than heat pump, kW
bodyThermalPower = gauss(0.1,0.5,K,t1); % thermal power from body heat, kW
solarThermalPower = 0.03*gauss(0.9,1.1,K,t1).*floorArea'.*solarFull; % thermal power from sunlight, kW 
qe = electricLoadThermalPower + bodyThermalPower + solarThermalPower; % exogenous thermal power, kW (8-12 W/m^2)

p1base = zeros(1,n1);
for k=1:K
    p1k = (((Tset(k+1,:) - a1.*Tbase2(k,:))./(1-a1) - thetaFull(k))./R + qe(k,:))./(0.78.*eta1(k,:)); % power to exactly track setpoint, kW
    p1k(p1k>0)=0;
    p1base(k,:) = min(pMax1(1,1),abs(p1k)); % saturated power, kW
    Tbase2(k+1,:) = a1.*Tbase2(k,:) + (1-a1).*(thetaFull(k) + R.*(eta1(k,:).*(-p1base(k,:)) + qe(k,:)));
end
     
electricPower = p1base;


waterHeatingLoad = zeros(K,n1);

if (~isnan(homeWithElectricWH) && homeWithElectricWH ~= 0)
waterHeatingLoad(:,1:round(n1*homeWithElectricWH)) = HPWHeachHouseLoad(:,1:round(n1*homeWithElectricWH));
end

% Calculate the cooling load values
coolingLoad = fullP + sum(electricPower,2) + sum(waterHeatingLoad,2) ;

summerPeak = quantile(coolingLoad,0.999);




end