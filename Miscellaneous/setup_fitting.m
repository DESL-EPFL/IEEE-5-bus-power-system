clear;

% data = load('simulation results');

t = data.t;
f = data.f(1,:);

t_cont = 100; % Time of contingency
idx_cont = find(t>=t_cont,1);

t = t(idx_cont:end);
f = f(idx_cont:end);

t = t-t(1);
df = f-f(1);