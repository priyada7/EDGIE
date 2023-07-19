function []=tranformerTempPlot(TOT1,TOT2,HST1,HST2,time,t,fs)
    figure(14)
    plot(time,TOT1,'LineWidth',3);
    xlabel('Hour','fontsize',fs)
    ylabel('Top oil temperature (C)','fontsize',fs)

    figure(15)
    plot(t(1:length(t)-1),HST1,'LineWidth',3);
    xlabel('Hour','fontsize',fs)
    ylabel('Hot spot temperature (C)','fontsize',fs)


   
    figure(16)
    plot(time,TOT2,'LineWidth',3);
    xlabel('Hour','fontsize',fs)
    ylabel('Top oil temperature (C)','fontsize',fs)

    figure(17)
    plot(t(1:length(t)-1),HST2,'LineWidth',3);
    xlabel('Hour','fontsize',fs)
    ylabel('Hot spot temperature (C)','fontsize',fs)
    end