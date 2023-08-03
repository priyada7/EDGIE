function []=tranformerTempPlot(TOT1,TOT2,HST1,HST2,time,t,fs)
% this function generates the temperature plots of the transformer oil
%
% Input:
%  TOT1, Kx1 vector for top-oil temperature unoptimized, C
%  TOT2, Kx1 vector for top-oil temperature unoptimized, C
%  HST1, Kx1 vector for hot-spot temperature unoptimized, C
%  HST2, Kx1 vector for hot-spot temperature optimized, C
%  time, timespan, h
%  t, (K+1)x1 vector for timespan, h
%  fs, fontsize
%
% Output:
%  top-oil temperature plots for optimized and unoptimized
%  hot-Spot temperature plots for optimized and unoptimized

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