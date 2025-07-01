function [Ad, Bd, Bd_d, n, mc, md] = dyn_heatpump()
% Given values
Cp_r = 810;   % kJ/°C
Cp_f = 3315;  % kJ/°C
Cp_w = 836;   % kJ/°C
UA_ra = 28;   % kJ/(°C h)
UA_fr = 624;  % kJ/(°C h)
UA_wf = 28;   % kJ/(°C h)
eta = 3;

% Coefficients
a11 = (-(UA_fr) - (UA_ra)) / Cp_r;
a12 = (UA_fr) / Cp_r;
a21 = (UA_fr) / Cp_f;
a22 = (-(UA_wf) - (UA_fr)) / Cp_f;
a23 = (UA_wf) / Cp_f;
a32 = (UA_wf) / Cp_w;
a33 = -(UA_wf) / Cp_w;

% Matrices
A = [a11, a12,  0;
     a21, a22, a23;
       0, a32, a33];

Bd = [UA_ra / Cp_r, 1 / Cp_r;
     0,           0;
     0,           0];

B = [0; 0; eta / Cp_w];

C = [1, 0, 0];


n = size(A,1); % Number of states
mc = size(B,2); % Number of controllable inputs
md = size(Bd,2); % Number of disturbance inputs

Ts = 0.25; %sampling time in hours (15 minutes)

% Augment system: Treat E as an additional input
sys_c = ss(A, [B Bd], C, 0);  % Combined input [B E]

% Discretize using Zero-Order Hold (ZOH)
sys_d = c2d(sys_c, Ts, 'zoh');

% Extract discrete-time matrices
Ad = sys_d.A;
Bd_Bdd = sys_d.B;  % [Bd Bdd] combined
Cd = sys_d.C;
Dd = sys_d.D;  % Should be zero

% Split Bd and Ed from the combined matrix
Bd = Bd_Bdd(:, 1);      % First column corresponds to Bd
Bd_d = Bd_Bdd(:, 2:end);  % Remaining columns correspond to Ed

end

