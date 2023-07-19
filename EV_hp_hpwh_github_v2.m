%% introduction
% This script simulates the grid impacts of electric vehicles, heat pumps and heat
% pump water heaters on the grid.

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% identifiers that users would like to simulate %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

batteryDeg =   1; % battery degradation with outside temperature is implemented
matpower   =   1; % set 1 if thermal model of transformer and voltage regulations needs to be solved
optimization = 1; % set 1 for v2h
rng(1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% input data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% timing (hour 0 is noon to make EV simulation easier)
ti = 0;                      % initial time, h
nDays=5;                     % no. of days
tf = nDays*24;               % final time, h
dt = 1;                      %  time step, h
K = tf/dt;                   % number of time steps
t = (0:dt:tf)';              % time span, h
n1 =33*2;                    % number of homes (= number of HPs)
carsPerHome = 2;             % number of cars per home
n2 = carsPerHome*n1;         % number of EVs
perBushome = 2;              % factor which decides how many house loads will be on each bus of 33 bus network (normally 990 homes gives 30 homes on each bus in a 33 bus network)
L =n1;                       % number of water heater

weatherfileName = 'cambridge-weather-2019.csv';                   % load weather file
loadfileName = 'MFRED-2019-NYC-Apartments-Electricity-Data.csv';  % load base load file
waterfile    = 'DHWEventGeneratorOutput.csv';                     % load water schduler file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% download weather & basline load data %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import weather data
weatherTime = (datetime(2019,1,1,0,0,0):hours(dt):datetime(2020,1,1,0,0,0))';
[thetaFull,solarFull] = importWeather(weatherfileName,weatherTime);
tStart = datetime(2019,1,18,0,0,0);           % Start date for simulation
tEnd = tStart + hours(tf);                    % End date for simulation
tt = (tStart:hours(dt):tEnd)';
theta = thetaFull(weatherTime>=tStart & weatherTime<tEnd);
solar = solarFull(weatherTime>=tStart & weatherTime<tEnd);

% import 'no electrification' load profile
[P,summerPeak] = importBaselineElectricity(loadfileName,tStart,tEnd,weatherTime,n1);
[pWorkBase] = workpower(n1,K,t);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% generate parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate heat pump parameters
[a1,C,R,eta1,Tset,qe,pMax1,pMaxAux]=generateBuildingModels(n1,dt,K,theta,t,P,solar);

% generate water heater parameters
[Rw,aw,ww,thetaw,Tsetw,etaw,pMaxHP,pMaxR] = generateWaterHeaterModels(waterfile,L,dt,tStart,tEnd);

% electric vehicle parameters
[a2,etac,etad,eMin,eMax,pcMax,pdMax,tau,e0,p0,ec,atHome,atWork,onRoad,tw,th] = generateEVmodels(n2,nDays,K,theta,batteryDeg,dt,t);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% baseline simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% heat pumps
[p1base,HPload,Tbase2]=heatPumpSimulation(K,n1,pMax1,eta1,pMaxAux,Tset,a1,R,qe,theta);

% electric vehicles (charging only at home, discharging only to drive)
[phBase,pwBase,eBase,pdBase]=evSimulation(K,n2,e0,p0,eMin,eMax,pcMax,a2,tau,etac,etad,dt,ec,onRoad,atHome,atWork);

 % heat-pump water heaters
[p1basehpwh,HPWHload,Tbasew,Tsetw]=heatPumpWaterHeaterSimulation(K,L,pMaxHP,etaw,pMaxR,aw,Rw,ww,thetaw);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% generate plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set scale of plot
lw = 2; % line width
fs = 24; % font size
tLim = [ti,tf]; % time axis limits, h
tTicks = ti:4:tf; % time axis ticks, h

if max(P)<1e3
    s = 1; % unit scalar, divide by s to convert from kW to plot units
    yMax = ceil(max(P/s + sum(phBase,2)/s+sum(p1base,2)/s)/500)*500; % y-axis limit
else
    s = 1e3; % for MW
    yMax = ceil(max(P/s + sum(phBase,2)/s+sum(p1base,2)/s)); % y-axis limit
end


baselinePlots(t,Tbase2,K,p1base,eBase,eMax,phBase,pwBase,pdBase,lw,fs,tLim,tTicks,tt,P,s,HPWHload,HPload,summerPeak)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if optimization ==1
    [e,ph,pw,pd,T,p1,p2,cvx_status]=optimizationCVX(eta1,pMax1,Tset,Tsetw,a1,K,theta,R,qe,pMaxAux,aw,ww,thetaw,pMaxR,a2,e0,eMax,eMin,tau,ec,dt,pdMax,pcMax,atWork,atHome,n1,n2,P,pWorkBase,Tbase2,Tbasew,Rw,etac,etad,t,tw,th);

    %optimized plots
    optimizedPlots(tt,t,T,K,p1,e,eMax,ph,pw,pd,lw,fs,tLim,tTicks,tf,P,p1base,phBase,atHome,atWork,summerPeak,pWorkBase,p1basehpwh,p2,s);

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% matpower & transformer modeling %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if matpower==1
    [IpuUnopt,IpuOpt,Tamb]=matpowerSimulation(phBase,ph,atHome,pd,n2,nDays,P,n1,p1basehpwh,p1base,p1,p2,perBushome,fs,theta);
    [TOT1, HST1, time] = transformer_git(Tamb, IpuUnopt);
    [TOT2, HST2, time] = transformer_git(Tamb, IpuOpt);
    tranformerTempPlot(TOT1,TOT2,HST1,HST2,time,t,fs);
    [lossUnoptT,lossOptT]=lossLife(HST1,HST2);
end
