    
function  [IpuUnopt,IpuOpt,Tamb]=matpowerSimulation(phBase,ph,atHome,pd,n2,nDays,...
    P,n1,p1basehpwh,p1base,p1,p2,perBushome,fs,theta)
 
% this function is used to generate compute power distribution calulation using MATPOWER 
% and generates the voltage profile for teh distribution network IEEE-33bus
%
% Input:
%  phBase, Kxn2 vector of power consumption from ev at home, kW
%  ph, Kxn2 matrix of home electrical consumption  due to ev, kW
%  atHome,indicator that vehicle's at home
%  pd, Kxn2 matrix of discharge power, kW
%  n2, number of ev
%  nDays, number of days
%  P, Kx1 vector of non-electrical load, kW
%  n1, number of home
%  p1basehpwh, Kx1 vector of water heater electrical load, kW
%  p1Base, Kxn1 vector of electrical consumption from heat pump, kW
%  p1, Kxn1 matrix of heat pump load, kW
%  p2, Kx1 vector of water heater electrical load, kW
%  perBushome, number of homes per bus
%  fs, fontsize parameter for plot
%  theta, Kx1 vector of outdoor temperature, C
% 
% Output:
%  IpuUnopt, Kx1 vector of unoptimized per unit current
%  IpuOpt, Kx1 vector of optimized per unit current
%  Tamb, Kx1 vector of outdoor temperature, C

    Unopt=phBase;
    Opt=ph-atHome.*pd;


    ind=1;
    BUnopt=zeros(nDays*24,n2);
    BOpt=zeros(nDays*24,n2);
    for i=1:(n2-1)
        BUnopt(:,i) =Unopt(:,ind)+Unopt(:,ind+1);
        BOpt(:,i)=Opt(:,ind)+Opt(:,ind+1);
        ind=ind+1;
    end

    j=(1:2:length(Unopt));
    CUnopt=zeros(nDays*24,n1);
    COpt=zeros(nDays*24,n1);
    for i=1:n1
        CUnopt(:,i) = BUnopt(:,j(i));
        COpt(:,i) = BOpt(:,j(i));
    end

    baseload =P'/n1;
    for i=1:n1
        baseload(i,:) = baseload(1,:);
    end
    loadprofileUnopt = (CUnopt + p1base + p1basehpwh(1:end-1,:))'+ baseload;
    loadprofileOpt = (COpt + p1 + p2 )'+ baseload;

    multiplier =1;
    matpowerloadUnopt= (multiplier*loadprofileUnopt);

    matpowerloadOpt= (multiplier*loadprofileOpt);



    beta = (1:perBushome:n1+perBushome);
    for i=1:length(beta)-1
        houseloadUnopt(i,:)= sum(matpowerloadUnopt(beta(:,i):beta(:,i+1)-1,:));
        houseloadOpt(i,:)= sum(matpowerloadOpt(beta(:,i):beta(:,i+1)-1,:));
    end

    mpcUnopt=loadcase('case33bw');
    mpcUnopt.gen(:,6)=1.04 ;  %substation voltage
  
    for i=1:120
        mpcUnopt.bus(:,3)=houseloadUnopt(:,i)/1e3;
        mpcUnopt.bus(:,4)=0.35.*houseloadUnopt(:,i)/1e3;
        resultsUnopt(:,i)=runpf(mpcUnopt);
        generationDataUnopt(:,i) = resultsUnopt(i).gen(:,2);
        loadDataUnopt = sum(houseloadUnopt)/1e3;
        LossDataUnopt(i,:) = generationDataUnopt(:,i) - loadDataUnopt(:,i);
        LossDataPercUnopt(i,:) = (generationDataUnopt(:,i)- loadDataUnopt(:,i))*100/generationDataUnopt(:,i);
    end

    mpcOpt=loadcase('case33bw');
    mpcOpt.gen(:,6)= mpcUnopt.gen(:,6);   %substation voltage
  
    for i=1:120
        mpcOpt.bus(:,3)=houseloadOpt(:,i)/1e3;
        mpcOpt.bus(:,4)=0.35.*houseloadOpt(:,i)/1e3;
        resultsOpt(:,i)=runpf(mpcOpt);
        generationDataOpt(:,i) = resultsOpt(i).gen(:,2);
        loadDataOpt = sum(houseloadOpt)/1e3;
        LossDataOpt(i,:) = generationDataOpt(:,i) - loadDataOpt(:,i);
        LossDataPercOpt(i,:) = (generationDataOpt(:,i)- loadDataOpt(:,i))*100/generationDataOpt(:,i);
    end



    figure(8)
    ylim([0.5 1.05]);
    [maxLoss, index1] = max(LossDataUnopt(:,1));
    plot(resultsUnopt(index1).bus(:,1),resultsUnopt(index1).bus(:,8));
    xlabel('Bus','fontsize',fs)
    ylabel('Voltage in p.u.','fontsize',fs)
   

    figure(9)
    ylim([0.5 1.05]);
    [maxLoss, index2] = max(LossDataOpt(:,1))
    plot(resultsOpt(index2).bus(:,1),resultsOpt(index2).bus(:,8));
    xlabel('Bus','fontsize',fs)
    ylabel('Voltage in p.u.','fontsize',fs)
   


    figure(10)
    plot(resultsUnopt(index1).bus(:,1),resultsUnopt(index1).bus(:,8),'LineWidth',2); hold on
    plot(resultsOpt(index2).bus(:,1),resultsOpt(index2).bus(:,8),'LineWidth',2);
    xlabel(" Bus ",'fontsize',fs);
    ylabel('Voltage in p.u.','fontsize',fs)
    legend("Unoptimized","Optimized")
    ylim([0 max(ylim)])
    hold off



    Irated = 36600;


    IactualUnopt = sum(houseloadUnopt)/0.240;
    IactualOpt = sum(houseloadOpt)/0.240;
    Tamb=theta;
    IpuUnopt = (round(IactualUnopt/Irated,1))';
    IpuOpt = (round(IactualOpt/Irated,1))';

end
