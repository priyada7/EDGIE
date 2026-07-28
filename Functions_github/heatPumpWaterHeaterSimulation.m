
function [HPWHload,HPWHeachHouseLoad]=heatPumpWaterHeaterSimulation(K,L,pMaxHP,etaw,pMaxR,aw,Rw,ww,thetaw,waterheater)

% this function performs the simulation of resistance water heater
%
% Input:
%  K, number of time steps
%  L, number of homes with water heaters
%  pMaxHP, 1xL vector of max electical capcity of heat pump, kW
%  etaw, KxL matrix of heat pump COP
%  pMaxR, 1xL vector of max resistance power, kW
%  aw, 1xL matrix of discrete-time dynamics parameters
%  Rw, 1xL vector of thermal resistances, C/kW
%  ww, KxL matrix of exogenous thermal power, kW
%  thetaw, Kx1 vector of water temperature, C
%
%
% Output:
%  p1basehpwh, KxL matrix for heat-pump water heater electric power, kW
%  HPWHload, Lx1 vector for heat-pump water heater electrical power, kW
%  Tbasew, (K+1)xL matrix for water temperature, C
%  Tsetw, (K+1)xL matrix of heating temperature setpoint, C

Kw=length(ww);
Qcw2=zeros(Kw,L);
Tsetw = repmat(f2c(trirnd(120,130,1,L)),Kw+1,1);
Tbasew = zeros(Kw+1,L); % water temperature, C
Tbasew(1,:) = Tsetw(1,:);
qw = zeros(Kw,L);      % heat pump electric power, kW
pw = zeros(Kw,L);          % resistance backup electric power, kW



p1basehpwh = zeros(Kw,L);
for k=1:Kw
 Qcw2 = (((Tsetw(k+1,:) - aw.*Tbasew(k,:))./(1-aw) - thetaw(k))./Rw + ww(k,:)); % power to exactly track setpoint, kW

if waterheater ==1
     
         qw(k,:) = max(0,min(pMaxR,Qcw2));
         pw(k,:) = qw(k,:);
     

 end


if waterheater ==2
     
         qw(k,:) = max(0,min(etaw(k,:).*pMaxHP(1,:),Qcw2));
         pw(k,:) = qw(k,:)./etaw(k,:);
     

 end

 if waterheater ==3
     if Qcw2 < etaw(k,:).*pMaxHP(1,:)
         qw(k,:) = max(0,min(etaw(k,:).*pMaxHP(1,:),Qcw2));
         pw(k,:) = qw(k,:)./etaw(k,:);
     else
         qw(k,:) = max(0,min(etaw(k,:).*pMaxHP(1,:) + pMaxR(1,:),Qcw2));
         pw(k,:) = max(qw(k,:)./etaw(k,:),qw(k,:)+(1-etaw(k,:)).*pMaxHP(1,:));
     end

 end
    Tbasew(k+1,:) = aw.*Tbasew(k,:) + (1-aw).*(thetaw(k) + Rw.*( qw(k,:) - ww(k,:) ));
    p1basehpwh = pw;
end
 HPWHeachHouseLoad = p1basehpwh;
 HPWHload = sum(p1basehpwh,2);

end
