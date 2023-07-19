function optimizedPlots(tt,t,T,K,p1,e,eMax,ph,pw,pd,lw,fs,tLim,tTicks,tf,P,p1base,phBase,atHome,atWork,summerPeak,pWorkBase,p1basehpwh,p2,s)

% temperature plot
figure(4), clf
subplot(2,1,1), plot(t,T,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks)
ylabel({'Indoor','temperature ($^\circ$C)'},'fontsize',fs)
title('Heat pumps','fontsize',fs)

% HP power plot
subplot(2,1,2), plot(t(1:K),p1,'linewidth',lw), grid on
hold on, plot(t(1:K),mean(p1,2),'k','linewidth',2*lw)
xlim(tLim), xticks(tTicks), ylim([0,max(ylim)])
ylabel({'Electric','power (kW)'},'fontsize',fs)
xlabel('Hour of day (0 = midnight)','fontsize',fs)

% individual energy
figure(5), clf
subplot(4,1,1), stairs(t,e./eMax,'linewidth',lw), grid on
hold on, stairs(t,mean(e./eMax,2),'k','linewidth',2*lw)
ylim([0,max(ylim)]), ylabel({'State of','charge'},'fontsize',fs)
xlim([0,tf]), xticks(tTicks)

% individual home charge power
subplot(4,1,2), plot(t(1:K),ph,'linewidth',lw), grid on
ylim([0,max(ylim)]), ylabel({'Home charge','power (kW)'},'fontsize',fs)
xlim([0,tf]), xticks(tTicks)

% individual work charge power
subplot(4,1,3), plot(t(1:K),pw,'linewidth',lw), grid on
ylim([0,max(ylim)]), ylabel({'Work charge','power (kW)'},'fontsize',fs)
xlim([0,tf]), xticks(tTicks)

% individual discharge power
subplot(4,1,4), plot(t(1:K),pd,'linewidth',lw), grid on
ylim([0,max(ylim)]), ylabel({'Discharge','power (kW)'},'fontsize',fs)
xlim([0,tf]), xticks(tTicks)
xlim([0,tf]), xticks(tTicks), xlabel('Hour','fontsize',fs)

figure(6), clf 
subplot(2,1,1),fullStairs(t,P + sum(p1base,2) + sum(phBase,2),'m','linewidth',lw), grid on
hold on, fullStairs(t,P + sum(p1,2) + sum(ph - atHome.*pd,2),'k','linewidth',lw)
stairs(t,0*t+summerPeak,'b--','linewidth',lw)
ylim([0,max(ylim)])
ylabel({'Aggregate home','power (kW)'},'fontsize',fs)

subplot(2,1,2), fullStairs(t,pWorkBase,'m','linewidth',lw), grid on
hold on, fullStairs(t,pWorkBase+sum(pw - atWork.*pd,2),'k','linewidth',lw)
ylim([0,max(ylim)])
ylabel({'Aggregate work','power (kW)'},'fontsize',fs)

figure(7), clf 
fullStairs(tt,P/s + sum(p1base,2)/s + sum(phBase,2)/s + sum(p1basehpwh(1:end-1,:),2)/s,'k','linewidth',lw), grid on
hold on, fullStairs(tt,P/s + sum(p1,2)/s + sum(ph - atHome.*pd,2)/s + sum(p2,2)/s ,'m--','linewidth',lw)
ylim([0,max(ylim)])
grid on
grid minor
if s==1e3, ylabel({'Aggregate',' power (MW)'},'fontsize',fs), end
if s==1, ylabel('Aggregate power (kW)','fontsize',fs), end
legend('Unoptimized','Optimized',Location="best")
end