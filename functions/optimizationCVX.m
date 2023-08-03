
function [e,ph,pw,pd,T,p1,p2,cvx_status]=optimizationCVX(eta1,pMax1,Tset,Tsetw,a1,K,theta,R,qe,pMaxAux,...
    aw,ww,thetaw,pMaxR,a2,e0,eMax,eMin,tau,ec,dt,pdMax,pcMax,atWork,atHome,...
    n1,n2,P,pWorkBase,Tbase2,Tbasew,Rw,etac,etad,t,tw,th)

% this function is used to compute optimization especially, vehicle-to-grid operation 
%
% Input:
%  eta1, Kxn1 matrix of heat pump COP
%  pMax1, kxn1 matrix of max electical capcity of heat pump, kW
%  Tset, (K+1)xn1 matrix of heating temperature setpoint, C
%  Tsetw, (K+1)xL matrix of heating temperature setpoint, C
%  a1, 1xn1 matrix of discrete-time dynamics parameters
%  K, number of time steps
%  theta, Kx1 vector of outdoor temperature, C
%  R, 1xn1 vector of thermal resistances, C/kW 
%  qe, (K+1)xn1 matrix of exogenous thermal power, kW 
%  pMaxAux, kxn1 matrix of max resistance power, kW
%  aw, 1xL matrix of discrete-time dynamics parameters
%  ww, KxL matrix of exogenous thermal power, kW
%  thetaw, Kx1 vector of water temperature, C
%  pMaxR, 1xL vector of max resistance power, kW
%  a2, 1xn2 vector of discrete-time dynamics parameters
%  e0, initial energy of ev battery, kWh
%  eMin, n2xK matrix of min-energy capacity, kWh 
%  eMax, n2xK matrix of max-energy capacity, kWh 
%  tau, KxL matrix of dissipation rate, 1/h
%  ec, 1xn2 vector of commute energy, kWh
%  dt, time step, h
%  pdMax, n2xK matrix of discharge capacity, kW
%  pcMax, n2xK matrix of charge capacity, kW
%  atHome,indicator that vehicle's at home
%  atWork, indicator that vehicle's at work
%  n1, number of homes
%  n2, number of ev
%  P, Kx1 vector of non-electrical load, kW
%  pWorkBase, Kx1 vector of woek place electrical load, kW
%  Tbase2, (K+1)xn1 matrix for indoor temperature, C
%  Tbasew, (K+1)xL matrix for water temperature, C
%  Rw, 1xL vector of thermal resistances, C/kW
%  etac, 1xn2 vector of charge efficiency
%  etad, 1xn2 vector of discharge efficiency
%  t,(K+1)x1 vector of timesteps, h
%  tw, 1xn2 vector for commute to work time of day, h
%  th, 1xn2 vector forcommute to home time of day, h
%
% Output:
%  e, (K+1)xn2 matrix of stored energy, kWh
%  ph, Kxn2 matrix of home electrical consumption  due to ev, kW
%  pw Kxn2 matrix of work electrical consumption  due to ev, kW
%  pd, Kxn2 matrix of discharge power, kW
%  T, (K+1)xn1 matrix of indoor temperature, h
%  p1, Kxn1 matrix of heat pump load, kW
%  p2, Kx1 vector of water heater electrical load, kW


pie = 0.15; % energy price, $/kWh
pid = 50; % peak demand price, $/kW
dT = 2; % allowable temperature deviation from setpoint, C
pic = 0.5; % discomfort price, $/C
%cvx_solver_settings('MIPGapAbs',1e-3,'MIPGap',1e-3,'NumericFocus',3)
cvx_begin quiet
     variables T(K+1,n1) Tw(K+1,n1) Qdotc(K,n1)  p2(K,n1)  e(K+1,n2) ph(K,n2) pw(K,n2) pd(K,n2)
     expressions p1(K,n1)  % intermediate variables (these are just placeholders for text expressions)
     p1 = Qdotc./eta1 + (1-1./eta1).*pos(Qdotc - eta1.*pMax1); % thermal -> electric power map (note assignment (=) vs. equality (==) operator)
   
     minimize( pid*max(P + sum(p1,2)+ sum(p2,2) + sum(ph - atHome.*pd,2))...        % home neighborhood peak
        + pid*max(0,max(pWorkBase + sum(pw - atWork.*pd,2)) - max(pWorkBase))...    % work neighborhood peak
        + pie*sum(P + sum(p1,2) + sum(p2,2) + sum(ph + pw,2))...                    % energy
        + pic*(sum(sum(abs(T-Tset))) -  sum(sum(abs(Tbase2-Tset)))) ...
        +  pic*(sum(sum(abs(Tw-Tsetw))) -  sum(sum(abs(Tbasew-Tsetw)))))            % discomfort
    % discomfort
 subject to                                                                         % constraints
%heat pump
 T(1,:) == Tset(1,:);                                                               % initial condition
 T(K+1,:) == Tset(K+1,:);                                                           % final condition: return to setpoint
 T(2:K+1,:) == repmat(a1,K,1).*T(1:K,:)...
     + (1-repmat(a1,K,1)).*(repmat(theta,1,n1) + repmat(R,K,1).*(Qdotc + qe));
 zeros(K,n1) <= Qdotc <= (eta1.*pMax1 + pMaxAux).*ones(K,n1)                        % heat pump electric power capacity limits
 abs(T - Tset) <= dT

 %heat-pump water heater
 Tw(1,:) == Tsetw(1,:);                                                             % initial condition
 Tw(K+1,:) == Tsetw(K+1,:);                                                         % final condition: return to setpoint
 Tw(2:K+1,:) == repmat(aw,K,1).*Tw(1:K,:)...
     + (1-repmat(aw,K,1)).*(repmat(thetaw,1,n1) + repmat(Rw,K,1).*(p2 - ww));
 abs(Tw - Tsetw) <= 10
 zeros(K,n1) <= p2 <= pMaxR.*ones(K,n1)


 % EV energy
    e(1,:) == e0                                                                    % initial condition
    e(K+1,:) == 0.3.*eMax(1,:)                                                      % terminal condition
    e(2:K+1,:) == repmat(a2,K,1).*e(1:K,:)...
        + (1-repmat(a2,K,1)).*repmat(tau,K,1).*...
        (repmat(etac,K,1).*(ph + pw) - pd./repmat(etad,K,1))
    repmat(eMin,K+1,1) <= e <= repmat(eMax,K+1,1)
    
    % EV power
    zeros(K,n2) <= ph <= repmat(pcMax,K,1)
    zeros(K,n2) <= pw <= repmat(pcMax,K,1)
    zeros(K,n2) <= pd <= repmat(pdMax,K,1)
    for i=1:n2
        pd(mod(t,24)==tw(i),i) == etad(i)*ec(i)/dt                                  % commute to work
        pd(mod(t,24)==th(i),i) == etad(i)*ec(i)/dt                                  % commute to home
       
    end
    ph(~atHome) == 0 % at-home charging
    pw(~atWork) == 0 % at-work charging
        
 
cvx_end
fprintf('CVX exit status: %s.\n',cvx_status)
% expand sparse matrices
e = full(e);
ph = full(ph);
pw = full(pw);
pd = full(pd);
T = full(T);
p1 = full(p1);
p2 = full(p2);


%%
if cvx_status == 'Failed'
    return
end
end