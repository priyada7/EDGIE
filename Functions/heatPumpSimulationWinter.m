function [p2,HPload,Tbase2,p_hp,aux]=heatPumpSimulationWinter(K,n1,cap2,eta2,pMaxAux,Tset,a1,R,qe,theta)
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

if n1==0

    p2 =0;
    HPload =0;
    Tbase2=0;
    p2 =0;
    p_hp =0;
    aux=0;

else
Tbase2 = zeros(K+1,n1);      % indoor temperature, C
Tbase2(1,:) = Tset(1,:);
q = zeros(K,n1);             % heat pump thermal power, kW
p2 = zeros(K,n1);             % electric power, kW
p_hp = zeros(K,n1); 
aux = zeros(K,n1); 
CapRentention = 0.8;

%cap2 =repmat(cap2,1,n1);



for k=1:K
    qck = (((Tset(k+1,:) - a1.*Tbase2(k,:))./(1-a1) - theta(k))./R - qe(k,:)); % power to exactly track setpoint, kW
    qck(qck<0)=0;
    %p2 =  max(q(theta)./cop(theta), cap(theta)./cop(theta) + qck - cap(theta));

    p2(k,:) =  max(qck./eta2(k,:), cap2(k,:)./eta2(k,:) + qck - cap2(k,:));
    p_hp(k,:) = min(qck ./ eta2(k, :), cap2(k, :) ./ eta2(k, :)); % Heat pump power (limited by capacity)
    aux(k,:) = max(0,qck - cap2(k,:));

    aux(aux>pMaxAux) = pMaxAux(1,1);

    idx = p2(k,:) > pMaxAux(k,:) + cap2(k,:) ./ eta2(k,:);
    p2(k,idx) = pMaxAux(k,idx) + cap2(k,idx) ./ eta2(k,idx);

    q(k,:) =  p_hp(k,:).* eta2(k, :) + aux(k,:);

    Tbase2(k+1,:) = a1.*Tbase2(k,:) + (1-a1).*(theta(k) + R.*(q(k,:)+ qe(k,:)));
end
     HPload = sum(p_hp+aux,2);
end
end