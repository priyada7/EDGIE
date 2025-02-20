function [p1base,HPload,Tbase2]=heatPumpSimulationWinter(K,n1,pMax1,eta2,pMaxAux,Tset,a1,R,qe,theta)
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

Tbase2 = zeros(K+1,n1);      % indoor temperature, C
Tbase2(1,:) = Tset(1,:);
p1base = zeros(K,n1);        % heat pump electric power, kW
q = zeros(K,n1);             % heat pump thermal power, kW
p = zeros(K,n1);             % electric power, kW

for k=1:K
    qck = (((Tset(k+1,:) - a1.*Tbase2(k,:))./(1-a1) - theta(k))./R - qe(k,:)); % power to exactly track setpoint, kW
    
    qck_condition = qck < eta2(k, :) .* pMax1(k, :);
   

% Case 1: qck_condition is true
q(k, qck_condition) = max(0, min(eta2(k, qck_condition) .* pMax1(k, qck_condition), qck(qck_condition)));
p(k, qck_condition) = q(k, qck_condition) ./ eta2(k, qck_condition);

% Case 2: qck_condition is false
q(k, ~qck_condition) = max(0, min(eta2(k, ~qck_condition) .* pMax1(k, ~qck_condition) + pMaxAux(k, ~qck_condition), qck(~qck_condition)));
p(k, ~qck_condition) = max(q(k, ~qck_condition) ./ eta2(k, ~qck_condition), q(k, ~qck_condition) + (1 - eta2(k, ~qck_condition)) .* pMax1(k, ~qck_condition));
    

    Tbase2(k+1,:) = a1.*Tbase2(k,:) + (1-a1).*(theta(k) + R.*(q(k,:)+ qe(k,:)));
    p1base = p;
end
     HPload = sum(p1base,2);
end