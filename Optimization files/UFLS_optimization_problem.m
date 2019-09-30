clear;
startup_UFLS;

%% Load simulation results
% load('Simulation results');

%% Extract data and parameters
h = UFLS.h;
N = UFLS.N;
sl = UFLS.sl;
pv = UFLS.pv;
pq = UFLS.pq;
Psp = UFLS.Psp;
Qsp = UFLS.Qsp;
Vsp = UFLS.Vsp;
fn = UFLS.fn;
Sn = UFLS.Sn;
Sgrid = UFLS.Sgrid;
Y = UFLS.Y;
YL = UFLS.YL;
YT = UFLS.YT;
kpv = UFLS.kpv;
kqv = UFLS.kqv;
kpf = UFLS.kpf;
kqf = UFLS.kqf;
Rp = UFLS.Rp(sort([sl;pv]));
H = UFLS.H;
D = UFLS.D;
T = UFLS.T;
fmin_tr = UFLS.fmin_tr/fn;
fmax_tr = UFLS.fmax_tr/fn;
fmin_st = UFLS.fmin_st/fn;
fmax_st = UFLS.fmax_st/fn;
Pmax = UFLS.Pmax;
Vmin = UFLS.Vmin;
Vmax = UFLS.Vmax;
Imax = UFLS.Imax; 
wrt_nodes = UFLS.wrt_nodes;
pv_no_sl = UFLS.pv_no_sl;
pq_no_sl = UFLS.pq_no_sl;
t_cont = UFLS.t_cont;
t_sim = UFLS.t_sim;

nsl = length(sl);
npv = length(pv);
npq = length(pq);
n = nsl + npv + npq;

% Tripped unit
tr = 2;
tr_no_sl = 1; % Index of tripped unit without slacks
LoG = 300/Sgrid; % Loss of generation (p.u.)

% Extract data at time of contingency
idx_cont = find(t >= t_cont,1);
f = f(1,idx_cont); % Consider frequency in first bus only, since it's the same everywhere in steady-state
V = V(:,idx_cont);
P = P(:,idx_cont);
Q = Q(:,idx_cont);

f0 = 1; % Nominal frequency in p.u.

%% Compute sensitivity coefficients
I = zeros(n,n);
for i = 1:n
    for j = 1:n
        I(i,j) = YL(i,j)*(V(i)-V(j))+YT(i,j)*V(i);
    end
end

Psp(pq) = P(pq).*(abs(V(pq))./Vsp(pq)).^(-kpv(pq));
Qsp(pq) = Q(pq).*(abs(V(pq))./Vsp(pq)).^(-kqv(pq));
[Kp,Kq,Hp,Hq] = sensCoeffs(Y,YL,YT,V,I,Psp,Qsp,Vsp,kpv,kqv,sl,pv,pq,wrt_nodes);

%% Declaration of symbolic decision variables
dPLS = sdpvar(npq,1,'full');
dQLS = sdpvar(npq,1,'full');
df = sdpvar(1,N+1,'full');
dr = sdpvar(1,N+1,'full');
dP = sdpvar(npv+npq,1,'full');
dQ = sdpvar(npv+npq,1,'full');
dPG_sl = sdpvar(nsl,1,'full');
dPG_pv = sdpvar(npv,1,'full');
dPG_sl_tmp = sdpvar(nsl,1,'full');
dPG_pv_tmp = sdpvar(npv,1,'full');
dPL = sdpvar(npq,1,'full');
dQG = sdpvar(npv,1,'full');
dQL = sdpvar(npq,1,'full');
dV = sdpvar(npv+npq,1,'full');
dI = sdpvar(n,n,'full');
delta = binvar(nsl+npv,1,'full'); % Binary decision variables for Droop characteristics
y = sdpvar(nsl+npv,N+1,'full'); % Contribution of generator to change of frequency

M1 = (fmax_tr-fmin_tr)./Rp.*Sn/Sgrid; % Big-M constants for piece-wise linear constraints related to variation of Droop constant
M2 = Pmax; % Big-M constants for piece-wise linear constraints related to Droop characteristics

%% Set initial values and compute equivalent system parameters
dr(1) = 0;
df(1) = 0;
y(:,1) = zeros(nsl+npv,1);

H_eq = sum(H.*Sn/Sgrid); % Equivalent system inertia
R_eq = 1/sum(Sn./(Rp*Sgrid)); % Equivalent system Droop

%% Objective function and constraints
objective = sum(dPLS);

cons = [];

% Cannot shed negative amount of load or more than 50% of actual load
cons = [cons, (zeros(npq,1) <= dPLS <= 0.5*abs(P(pq)))];

% Estimation of shed reactive power (assuming constant cos(phi)
cons = [cons, (dQLS == dPLS.*Q(pq)./P(pq))];

%% Run through the trajectory
for i = 1:N
    %% Frequency variation
    cons = [cons, (-delta.*M1 + Sn./Rp*df(i+1)/Sgrid <= y(:,i+1) <= Sn./Rp*df(i+1)/Sgrid + delta.*M1)];
    cons = [cons, (-(1-delta).*M1 <= y(:,i+1) <= (1-delta).*M1)];
    cons = [cons, (df(i+1) == df(i) + h*f0/(2*H_eq)*(dr(i+1)-LoG+sum(dPLS)-D*df(i+1)))];
    cons = [cons, (dr(i+1) == dr(i) + h*1/T*(-dr(i+1)-sum(y(:,i+1))))];

    cons = [cons, (fmin_tr <= f + df(i+1) <= fmax_tr)];
end

%% Compound vectors of generators and loads powers variations
for i = 1:npv
    if pv(i) == tr % Substract loss of generation in tripped unit
        cons = [cons, (dP(tr_no_sl) == dPG_pv(i)-LoG)];
    else
        cons = [cons, (dP(pv_no_sl(i)) == dPG_pv(i))];
    end
end

cons = [cons, (dQ(pv_no_sl) == dQG)];
cons = [cons, (dP(pq_no_sl) == dPL)];
cons = [cons, (dQ(pq_no_sl) == dQL)];

%% Gen power variation
for k = 1:nsl
    dPG_sl_tmp(k) = -1/Rp(sl(k))*df(end);
    cons = [cons, (dPG_sl_tmp(k) <= Pmax(sl(k))-P(sl(k))+delta(sl(k))*M2(sl(k)))];
    cons = [cons, (dPG_sl_tmp(k) >= Pmax(sl(k))-P(sl(k))-(1-delta(sl(k)))*M2(sl(k)))];
    cons = [cons, (-delta(sl(k))*M2(sl(k))+dPG_sl_tmp(k) <= dPG_sl(k) <= dPG_sl_tmp(k)+delta(sl(k))*M2(sl(k)))];
    cons = [cons, (-(1-delta(sl(k)))*M2(sl(k))+Pmax(sl(k))-P(sl(k)) <= dPG_sl(k) <= Pmax(sl(k))-P(sl(k))+(1-delta(sl(k)))*M2(sl(k)))];
end

for k = 1:npv
    dPG_pv_tmp(k) = -1/Rp(pv(k))*df(end);
    cons = [cons, (dPG_pv_tmp(k) <= Pmax(pv(k))-P(pv(k))+delta(pv(k))*M2(sl(k)))];
    cons = [cons, (dPG_pv_tmp(k) >= Pmax(pv(k))-P(pv(k))-(1-delta(pv(k)))*M2(sl(k)))];
    cons = [cons, (-delta(pv(k))*M2(sl(k))+dPG_pv_tmp(k) <= dPG_pv(k) <= dPG_pv_tmp(k)+delta(pv(k))*M2(sl(k)))];
    cons = [cons, (-(1-delta(pv(k)))*M2(sl(k))+Pmax(pv(k))-P(pv(k)) <= dPG_pv(k) <= Pmax(pv(k))-P(pv(k))+(1-delta(pv(k)))*M2(sl(k)))];
end

% Gen reactive power variation (pv generators only)
cons = [cons, (dQG == dP(pv_no_sl).*Q(pv)./P(pv))]; 

%% Load power variation
cons = [cons, (dPL == dPLS + Psp(pq).*kpv(pq)./Vsp(pq).*dV(pq_no_sl) + Psp(pq).*kpf(pq)*df(end))];
cons = [cons, (dQL == dQLS + Qsp(pq).*kqv(pq)./Vsp(pq).*dV(pq_no_sl) + Qsp(pq).*kqf(pq)*df(end))];   
                                
%% Voltage and current variation
cons = [cons, (dV == Kp*dP + Kq*dQ)];

dI_tmp = zeros(n,n);
for j = 1:length(wrt_nodes)
    dI_tmp = dI_tmp + Hp(:,:,j)*dP(j) + Hq(:,:,j)*dQ(j); 
end

cons = [cons, (dI == dI_tmp)];

%% Safety constraints
cons = [cons, (Vmin(sort([pv;pq])) <= abs(V(sort([pv;pq]))) + dV <= Vmax(sort([pv;pq])))]; 
cons = [cons, (abs(I) + dI <= Imax)];

% Steady-state frequency
cons = [cons, (fmin_st <= f + df(end) <= fmax_st)];

%% Solve the problem
options = sdpsettings('verbose',1,'debug',1,'solver','gurobi');
sol = optimize(cons,objective,options);

dPLS = value(dPLS);
dQLS = value(dQLS);          
df = value(df);
dP_tmp = value(dP);
dP = zeros(n,1);
dP(sl) = value(dPG_sl); % Include variation of active power in slacks
dP(pv) = dP_tmp(pv_no_sl);
dP(pq) = dP_tmp(pq_no_sl);
dQ = value(dQ);
dV = value(dV);
dI = value(dI);

%% Save results
clearvars -except UFLS dPLS dQLS df dP dQ dV dI
% save optimization results
