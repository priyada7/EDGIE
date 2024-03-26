
function [p1basehpwh,HPWHload,Tbasew,Tsetw]=heatPumpWaterHeaterSimulation(K,L,pMaxHP,etaw,pMaxR,aw,Rw,ww,thetaw,waterheater)

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

Kw=K;
Qcw2=zeros(Kw,L);
Tsetw = repmat(f2c(trirnd(120,130,1,L)),Kw+1,1);
Tbasew = zeros(Kw+1,L); % water temperature, C
Tbasew(1,:) = Tsetw(1,:);
qw = zeros(K,L);      % heat pump electric power, kW
pw = zeros(K,L);          % resistance backup electric power, kW


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
     Qcw2_condition = Qcw2 < etaw(k, :) .* pMaxHP(1, :);


% Case 1: Qcw2_condition is true
qw(k, Qcw2_condition) = max(0, min(etaw(k, Qcw2_condition) .* pMaxHP(1, Qcw2_condition), Qcw2(Qcw2_condition)));
pw(k, Qcw2_condition) = qw(k, Qcw2_condition) ./ etaw(k, Qcw2_condition);

% Case 2: Qcw2_condition is false
qw(k, ~Qcw2_condition) = max(0, min(etaw(k, ~Qcw2_condition) .* pMaxHP(1, ~Qcw2_condition) + pMaxR(1, ~Qcw2_condition), Qcw2(~Qcw2_condition)));
pw(k, ~Qcw2_condition) = max(qw(k, ~Qcw2_condition) ./ etaw(k, ~Qcw2_condition), qw(k, ~Qcw2_condition) + (1 - etaw(k, ~Qcw2_condition)) .* pMaxHP(1, ~Qcw2_condition));

 end
    Tbasew(k+1,:) = aw.*Tbasew(k,:) + (1-aw).*(thetaw(k) + Rw.*( qw(k,:) - ww(k,:) ));
    p1basehpwh = pw;
end

 HPWHload = sum(p1basehpwh,2);

end
% 
%     [row1, col1]=size(p1basew);
%     for i = 1 : row1
%         Qcw2(i,:) = ((Tsetw(i+1,:) - aw.*Tbasew(i,:))./(1-aw) - thetaw(i))./Rw + ww(i,:); % power to exactly track setpoint, kW
%         for j = 1 : col1
%             if Qcw2(i,j) <=0
%                 Qcw2(i,j)=0;
%                 pAuxhpwh(i,j) = 0;
%             elseif Qcw2(i,j) <= etaw(i,j)*pMaxHP(1,j)
%                 pAuxhpwh(i,j) = Qcw2(i,j);
%             elseif Qcw2(i,j) >etaw(i,j)*pMaxHP(1,j) & Qcw2(i,j)-etaw(i,j)*pMaxHP(1,j)<= pMaxR(1,j)
%                 pAuxhpwh(i,j) = Qcw2(i,j)/1;
%             else
%                 pAuxhpwh(i,j) = pMaxR(1,j);
%             end
%         end
%     end
%   
% 
%   Tbasew(k+1,:) = aw.*Tbasew(k,:) + (1-aw).*(thetaw(k) + Rw.*( pAuxhpwh(k,:) - ww(k,:) ));
%  p1basehpwh = pAuxhpwh;
% end

%HPWHload = sum(pAuxhpwh,2);
%end