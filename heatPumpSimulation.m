function [p1base,HPload,Tbase2]=heatPumpSimulation(K,n1,pMax1,eta1,pMaxAux,Tset,a1,R,qe,theta)

Tbase2 = zeros(K+1,n1);      % indoor temperature, C
Tbase2(1,:) = Tset(1,:);
p1base = zeros(K,n1);        % heat pump electric power, kW
pHpbase2 = zeros(K,n1);      % heat pump electric power, kW
pAux = zeros(K,n1);          % resistance backup electric power, kW

for k=1:K  
%  % to compute power input from Heat Pump
[row1, col1]=size(pHpbase2);
for i = 1 : row1
      Qck(i,:) = (((Tset(i+1,:) - a1.*Tbase2(i,:))./(1-a1) - theta(i))./R - qe(i,:)); 
      for j = 1 : col1
          if Qck(i,j) <= 0
             Qhp(i,j) =0;
             Qck(i,j) =0;
         elseif Qck(i,j) <= eta1(i,j)*pMax1(i,j) & Qck(i,j) >=0
          Qhp(i,j) = Qck(i,j);
          else
              Qhp(i,j) = eta1(i,j)*pMax1(i,j);
          end
          pHpbase2(i,j) = Qhp(i,j)/eta1(i,j);
      end
end
 % to compute power input from Aux
for i = 1 : row1
     Qck2(i,:) = (((Tset(i+1,:) - a1.*Tbase2(i,:))./(1-a1) - theta(i))./R - qe(i,:)); 
     for j = 1 : col1
          if Qck2(i,j) <=0
              Qck2(i,j) =0;
               pAux(i,j) = 0;
          elseif Qck2(i,j) <= eta1(i,j)*pMax1(i,j)
              pAux(i,j) = 0; 
         elseif Qck2(i,j)> eta1(i,j)*pMax1(i,j) & Qck2(i,j)- eta1(i,j)*pMax1(i,j) <= pMaxAux(i,j)
              pAux(i,j) = Qck2(i,j)-eta1(i,j)*pMax1(i,j);
          else
              pAux(i,j) = pMaxAux(i,j);
      end

      end
end
  Tbase2(k+1,:) = a1.*Tbase2(k,:) + (1-a1).*(theta(k) + R.*(eta1(k,:).*pHpbase2(k,:) +  pAux(k,:) + qe(k,:)));
   p1base = pHpbase2+pAux;
end
HPload   = sum(pHpbase2,2) + sum(pAux,2);


end