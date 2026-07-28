function [electricPowerWinter,waterHeatingLoad]=winterPeakcalc(zone,thetaFull,solarFull,dt,fullP,n1,weatherTime,R,cap2,HPWHloadFullyear,homeWithHP,HPWHeachHouseLoad,floorArea,pMaxAux,homeWithAux,homeWithElectricWH,percentageAt_DT,check)
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


homeWithAux = round(n1*percentageAt_DT*homeWithAux);
homeWithHP  = round(n1*percentageAt_DT*homeWithHP);
K=length(thetaFull);




if homeWithHP==0
    p1baseHP=0;
else
    Rhp=R(1:homeWithHP);
    K=length(thetaFull);
    Tbase2 = zeros(K+1,homeWithHP);      % indoor temperature, C
    Tset = repmat(f2c(trirnd(65,74,1,homeWithHP)),K+1,1);
    Tbase2(1,:) = Tset(1,:);

    floorAreahp = floorArea(1:homeWithHP);                   % conditioned floor area, m^2
    ceilingHeight = repmat(3,1,homeWithHP);
    rhoc = 3.42e-4;                                         % volumetric thermal capacitance of air, kWh/C/m^3
    C = trirnd(16,18,1,homeWithHP).*(rhoc*floorAreahp'.*ceilingHeight); % thermal capacitance, kWh/C
    a1 = exp(-dt./(Rhp.*C));

    thetaLow = f2c(5);                                      % first temperature point for heat pump COP curve, C
    thetaHigh = f2c(45); 

    if zone > 4
        etaLow = trirnd(1.5,1.8,1,homeWithHP);                            % first COP point
        etaHigh = etaLow + 1.5 + trirnd(0,0.5,1,homeWithHP);            % second COP point
        eta1 = etaLow + (etaHigh-etaLow).*(thetaFull-thetaLow)./(thetaHigh-thetaLow); % COP curve
        eta1(eta1 < 1)=1;
    else
        etaLow = trirnd(1.1,1.5,1,homeWithHP);                            % first COP point
        etaHigh = etaLow + 1.5 + trirnd(0,0.5,1,homeWithHP);            % second COP point
        eta1 = etaLow + (etaHigh-etaLow).*(thetaFull-thetaLow)./(thetaHigh-thetaLow); % COP curve
        eta1(eta1 < 1)=1;
    end

    cap2 = cap2(:,1:homeWithHP);

    electricLoadThermalPower = fullP/n1; % thermal power from electrical loads other than heat pump, kW
    bodyThermalPower = trirnd(0.1,0.5,K,homeWithHP); % thermal power from body heat, kW
    solarThermalPower = 0.03*trirnd(0.9,1.1,K,homeWithHP).*floorAreahp'.*solarFull; % thermal power from sunlight, kW
    qe = electricLoadThermalPower + bodyThermalPower + solarThermalPower; % exogenous thermal power, kW (8-12 W/m^2)
if check ==1
    qe=qe./4;
end

    q = zeros(K,homeWithHP);             % heat pump thermal power, kW
    p = zeros(K,homeWithHP);             % electric power, kW


    p1baseHP = zeros(1,homeWithHP);
    p1baseAux = zeros(1,homeWithAux);
    for k=1:K
        qck = (((Tset(k+1,:) - a1.*Tbase2(k,:))./(1-a1) - thetaFull(k))./Rhp - qe(k,:)); % power to exactly track setpoint, kW
   p(k,:) = min(qck ./ eta1(k, :), cap2(k, :) ./ eta1(k, :)); % Heat pump power (limited by capacity)
   q(k,:) =  p(k,:).* eta1(k, :);
   Tbase2(k+1,:) = a1.*Tbase2(k,:) + (1-a1).*(thetaFull(k) + Rhp.*(q(k,:)+ qe(k,:)));
   p1baseHP = p;

    end
end
%%%%%%

if homeWithAux==0
    p1baseAux=0;
else
    q = zeros(K,homeWithAux);             % heat pump thermal power, kW
    p = zeros(K,homeWithAux);             % electric power,

    Tbase2aux = zeros(K+1,homeWithAux);      % indoor temperature, C
    Tsetaux = repmat(f2c(trirnd(65,74,1,homeWithAux)),K+1,1);
    Tbase2aux(1,:) = Tsetaux(1,:);
   
    if homeWithHP ~= 0
        Raux=R(homeWithHP:homeWithHP+length(homeWithAux)-1);
        floorAreaaux = floorArea(homeWithHP:homeWithHP+length(homeWithAux)-1);                   % conditioned floor area, m^2
    else
        Raux=R(homeWithHP+1:homeWithHP+length(homeWithAux));
        floorAreaaux = floorArea(homeWithHP+1:homeWithHP+length(homeWithAux));
    end
    ceilingHeightaux = repmat(3,1,homeWithAux);
    rhoc = 3.42e-4;                                         % volumetric thermal capacitance of air, kWh/C/m^3
    Caux = trirnd(16,18,1,homeWithAux).*(rhoc*floorAreaaux'.*ceilingHeightaux); % thermal capacitance, kWh/C
    a1aux = exp(-dt./(Raux.*Caux));


    electricLoadThermalPower = fullP/n1; % thermal power from electrical loads other than heat pump, kW
    bodyThermalPower = trirnd(0.1,0.5,K,homeWithAux); % thermal power from body heat, kW
    solarThermalPower = 0.03*trirnd(0.9,1.1,K,homeWithAux).*floorAreaaux'.*solarFull; % thermal power from sunlight, kW
    qeAux = electricLoadThermalPower + bodyThermalPower + solarThermalPower; % exogenous thermal power, kW (8-12 W/m^2)

if check ==1
    qeAux=qeAux./4;
end
    pMaxAux = pMaxAux(:,1:homeWithAux);

    for k=1:K
        qck = (((Tsetaux(k+1,:) - a1aux.*Tbase2aux(k,:))./(1-a1aux) - thetaFull(k))./Raux - qeAux(k,:)); % power to exactly track setpoint, kW
        q(k,:) = max(0, min(pMaxAux(k,:), qck));
        p(k, :) = q(k, :) ./ 1;
        Tbase2aux(k+1,:) = a1aux.*Tbase2aux(k,:) + (1-a1aux).*(thetaFull(k) + Raux.*(q(k,:)+ qeAux(k,:)));
        p1baseAux = p;


  

    end
end
electricPowerWinter = sum(p1baseHP,2)+sum(p1baseAux,2);

waterHeatingLoad = zeros(K,n1);
if (~isnan(homeWithElectricWH) && homeWithElectricWH ~= 0)
    waterHeatingLoad(:,1:round(n1*homeWithElectricWH)) = HPWHeachHouseLoad(:,1:round(n1*homeWithElectricWH));
end
end