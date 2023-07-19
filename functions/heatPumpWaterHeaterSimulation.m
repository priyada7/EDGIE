
function [p1basehpwh,HPWHload,Tbasew,Tsetw]=heatPumpWaterHeaterSimulation(K,L,pMaxHP,etaw,pMaxR,aw,Rw,ww,thetaw)
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