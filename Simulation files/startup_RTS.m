clear;

%% Base values
fn = 50;
Vn = 345e3;
Sn = 1000e6;
Zn = Vn^2/Sn;
In = Sn/(sqrt(3)*Vn);

t_cont = 100; % Time of contingency

%% Total r, l, b in pu
r12 = 0.02;
x12 = 0.06;
b12 = 0.06;

r13 = 0.08;
x13 = 0.24;
b13 = 0.05;

r23 = 0.06;
x23 = 0.18;
b23 = 0.04;

r24 = 0.06;
x24 = 0.18;
b24 = 0.04;

r25 = 0.04;
x25 = 0.12;
b25 = 0.03;

r34 = 0.01;
x34 = 0.03;
b34 = 0.02;

r45 = 0.08;
x45 = 0.24;
b45 = 0.05;

% Order: 1    2    3    4    5    6    7
%       1-2, 1-3, 2-3, 2-4, 2-5, 3-4, 4-5
r = [r12;r13;r23;r24;r25;r34;r45];
x = [x12;x13;x23;x24;x25;x34;x45];
b = [b12;b13;b23;b24;b25;b34;b45];

%% Total R, L, C in Ohm, H, F
R = r*Zn;
L = x*Zn/(2*pi*60);
C = b/(2*pi*60*Zn);

%% Compute line lengths in order to achieve a propagation speed around 95% of light speed
T_prop = sqrt(L.*C);
len = ceil(0.95*300000*T_prop);

%% Zero-sequence parameters
R0 = 3*R;
L0 = 3*L;
C0 = C/3;

%% Load power profiles (in W)
PL1 = 300e6;
QL1 = 30e6;
PL2 = 300e6;
QL2 = 30e6;
PL3 = 300e6;
QL3 = 30e6;

%% Load shedding (in W)
dPLS1 = 0;
dQLS1 = 0;
dPLS2 = 0;
dQLS2 = 0;
dPLS3 = 0;
dQLS3 = 0;