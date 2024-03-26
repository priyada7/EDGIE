function []=baselinePlots(t,Tbase2,K,p1base,eBase,eMax,phBase,pwBase,pdBase,lw,fs,tLim,tTicks,tt,P,s,HPWHload,HPload,summerPeak)

% this function is used to generate base line plots and generates the
% overall load profile with different appliances 
%
% Input:
%  t,(K+1)x1 vector of timesteps, h
%  Tbase2, (K+1)xn1 vector of indoor temperature, C
%  p1Base, Kxn1 vector of electrical consumption from heat pump, kW
%  eBase, (K+1)xn2 vector of stored energy in battery, kWh
%  eMax, n2xK vector of energy capacity, kWh 
%  phBase, Kxn2 vector of power consumption from ev at home, kW
%  pwBase, Kxn2 vector of power consumption from ev at work, kW
%  pdBase, Kxn2 vector of ev discharge power, kW
%  lw, linedwidth parameter for plot
%  fs, fontsize parameter for plot
%  tLim, 1x2 vector of time limits
%  tt, date timespan
%  P, Kx1 vector of non-electrical load, kW
%  s, parameter to determine scale of plots
%  HPWHload, (K+1)xn1 vector of water heater electrical consumption, kW
%  HPload, Kxn1 vector of heatpump electrical consumption, kW
%  summerPeak, paramter for summer peak
%
% Output:
%  emperature plot
%  heat pump plot
%  energy plot
%  home charge plot
%  work charge plot
%  discharge plot
%  overall loadprofile plot

% temperature plot
figure(1), clf
subplot(2,1,1), plot(t,Tbase2,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks)
ylabel({'Indoor','temperature ($^\circ$C)'},'fontsize',fs)
title('Heat pumps','fontsize',fs)

% HP power plot
subplot(2,1,2), plot(t(1:K),p1base,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks), ylim([0,max(ylim)])
ylabel({'Electric','power (kW)'},'fontsize',fs)
xlabel('Hour','fontsize',fs)

% energy plot
figure(2), clf
subplot(4,1,1), plot(t,eBase./eMax,'linewidth',lw), grid on
hold on, plot(t,mean(eBase./eMax,2),'k','linewidth',2*lw)
xlim(tLim), xticks(tTicks)
ylabel({'State of','charge'},'fontsize',fs), ylim([0,1])
title('Electric vehicles','fontsize',fs)

% home charge plot
subplot(4,1,2), plot(t(1:K),phBase,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks), ylim([0,max(ylim)])
ylabel({'Home charge','power (kW)'},'fontsize',fs)

% work charge plot
subplot(4,1,3), plot(t(1:K),pwBase,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks), ylim([0,max(ylim)])
ylabel({'Work charge','power (kW)'},'fontsize',fs)

% discharge plot
subplot(4,1,4), plot(t(1:K),pdBase,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks), ylim([0,max(ylim)])
ylabel({'Discharge','power (kW)'},'fontsize',fs)
xlabel('Hour','fontsize',fs)

figure(3), clf
fullStairs(tt,P/s,LineWidth=2); hold on
fullStairs(tt,P/s + HPWHload/s,'m',LineWidth=2),
fullStairs(tt,P/s + HPWHload/s + sum(phBase,2)/s ,'g',LineWidth=2),
fullStairs(tt,P/s + HPWHload/s + sum(phBase,2)/s + HPload/s,'k',LineWidth=2)
yline(summerPeak/s,'r',LineWidth=2,LineStyle=":")
if s==1e3, ylabel('Aggregate power (MW)','fontsize',24), end
if s==1, ylabel('Aggregate power (kW)','fontsize',24), end
ylim([0 max(ylim)])
grid on
grid minor
legend('Base Electrified Load','Add HPWH ','Add EV','Add HP','Summer Peak',Location='best')

end