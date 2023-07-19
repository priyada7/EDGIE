

function [phBase,pwBase,eBase,pdBase]=evSimulation(K,n2,e0,p0,eMin,eMax,pcMax,a2,tau,etac,etad,dt,ec,onRoad,atHome,atWork)

eBase = zeros(K+1,n2); % stored energy, kWh
eBase(1,:) = e0;
pcBase = zeros(K,n2); % home or work charge power, kW
pdBase = zeros(K,n2); % discharge power, kW
for k=1:K
    % initialize charging power with previous chargin power
    if k==1
        pck = p0;
    else
        pck = pcBase(k-1,:);
    end
    
    % start charging if energy dropped below minimum
    pck(eBase(k,:)<eMin) = pcMax(eBase(k,:)<eMin);
    
    % stop charging if driving
    pck(~atHome(k,:) & ~atWork(k,:)) = 0;
    
    % saturate charging power
    pcBase(k,:) = max(0,min(pck,(eMax-a2.*eBase(k,:))./(1-a2)./tau./etac));
    
    % discharge power if on the road
    pdBase(k,onRoad(k,:)) = etad(onRoad(k,:)).*ec(onRoad(k,:))/dt;
    
    % update stored energy
    eBase(k+1,:) = a2.*eBase(k,:) + (1-a2).*tau.*(etac.*pcBase(k,:) - pdBase(k,:)./etad);
end
phBase = atHome.*pcBase; % baseline home charge power, kW
pwBase = atWork.*pcBase; % baseline work charge power, kW



end