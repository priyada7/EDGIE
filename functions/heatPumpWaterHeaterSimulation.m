
function [p1basehpwh,HPWHload,Tbasew,Tsetw]=heatPumpWaterHeaterSimulation(K,L,pMaxHP,etaw,pMaxR,aw,Rw,ww,thetaw)

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


pAuxhpwh= zeros(Kw+1,L);
p1basew = zeros(Kw,L); % electric power, kW
p1basehpwh = zeros(Kw,L);
for k=1:Kw
    [row1, col1]=size(p1basew);
    for i = 1 : row1
        Qcw2(i,:) = ((Tsetw(i+1,:) - aw.*Tbasew(i,:))./(1-aw) - thetaw(i))./Rw + ww(i,:); % power to exactly track setpoint, kW
        for j = 1 : col1
            if Qcw2(i,j) <=0
                Qcw2(i,j)=0;
                pAuxhpwh(i,j) = 0;
            elseif Qcw2(i,j) <= etaw(i,j)*pMaxHP(1,j)
                pAuxhpwh(i,j) = Qcw2(i,j);
            elseif Qcw2(i,j) >etaw(i,j)*pMaxHP(1,j) & Qcw2(i,j)-etaw(i,j)*pMaxHP(1,j)<= pMaxR(1,j)
                pAuxhpwh(i,j) = Qcw2(i,j)/1;
            else
                pAuxhpwh(i,j) = pMaxR(1,j);
            end
        end
    end
  

  Tbasew(k+1,:) = aw.*Tbasew(k,:) + (1-aw).*(thetaw(k) + Rw.*( pAuxhpwh(k,:) - ww(k,:) ));
 p1basehpwh = pAuxhpwh;
end

HPWHload = sum(pAuxhpwh,2);


end