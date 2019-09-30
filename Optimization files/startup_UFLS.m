clear;

% Add paths to working directory for Yalmip and Gurobi
% addpath(genpath('Path to Yalmip installation'));
% addpath(genpath('Path to Gurobi installation'));

% Load linedata and busdata files
linedata = load('linedata.txt');
busdata = load('busdata.txt');

[sl,pv,pq,Psp,Qsp,Vsp,Pmax,~,~,Vmin,Vmax,Imax,Y,YL,YT,kpv,kqv,kpf,kqf] = createGridFromFiles(busdata,linedata);

nsl = length(sl);
npv = length(pv);
npq = length(pq);
n = nsl+npv+npq;
n_lines = size(linedata,1);

fn = 50;

PG_loss = 0.025*Pmax; % Losses in generating unit (2.5% of nominal power)

% Prediction horizon and number of steps
t_sim = 200; % sec
h = 0.1; % Step
N = t_sim/h; % Number of steps

% Create UFLS structure for UFLS optimization problem
UFLS.h = h;
UFLS.N = N;
UFLS.sl = sl;
UFLS.pv = pv;
UFLS.pq = pq;
UFLS.Psp = Psp;
UFLS.Qsp = Qsp;
UFLS.Vsp = Vsp;
UFLS.fn = fn;
UFLS.Sn = [1000;1000]; % MVA. Nominal powers of the generators
UFLS.Sgrid = 1000; % MVA. Nominal power of the grid
UFLS.n_lines = n_lines;
UFLS.Y = Y;
UFLS.YL = YL;
UFLS.YT = YT;
UFLS.kpv = kpv;
UFLS.kqv = kqv;
UFLS.kpf = kpf;
UFLS.kqf = kqf;
UFLS.Rp = [0.05;0.05];
UFLS.H = [3.48;3.48]; % sec
UFLS.D = 2.67;
UFLS.T = 205; % sec
UFLS.fmin_tr = 47.5;
UFLS.fmax_tr = 52;
UFLS.fmin_st = 49.5;
UFLS.fmax_st = 50.5;
UFLS.Pmax = Pmax-PG_loss;
UFLS.Vmin = Vmin;
UFLS.Vmax = Vmax;
UFLS.Imax = Imax;
UFLS.wrt_nodes = sort([pv;pq]);
UFLS.t_cont = 100; % Time of contingency
UFLS.t_sim = t_sim; % Prediction horizon

% Remove slacks in pv, pq indices
pq_no_sl = pq;
pv_no_sl = pv;
for j = 1:n
    sl_idx = find(sl == j);
    if ~isempty(sl_idx) % i is slack. Following pv, pq indices have to be decremented  
        pq_tmp = find(pq > j);
        pv_tmp = find(pv > j);
        if ~isempty(pq_tmp)
            pq_no_sl(pq_tmp(:)) = pq_no_sl(pq_tmp(:)) - ones(length(pq_tmp),1);
        end
        if ~isempty(pv_tmp)
            pv_no_sl(pv_tmp(:)) = pv_no_sl(pv_tmp(:)) - ones(length(pv_tmp),1);
        end
    end
end

UFLS.pv_no_sl = pv_no_sl;
UFLS.pq_no_sl = pq_no_sl;

clearvars -except UFLS