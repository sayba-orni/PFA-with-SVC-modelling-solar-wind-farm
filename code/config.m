function C = config()
C.Sbase   = 100;
C.Vband   = [0.95 1.05];
C.Vref    = 1.00;

C.SVC_Q_MAX = 300;                          % MVAr rating
C.Bmax =  C.SVC_Q_MAX/(C.Vref^2*C.Sbase);
C.Bmin = -C.Bmax;

% 5-bus network
C.pairs = [1 2; 1 3; 2 3; 2 4; 2 5; 3 4; 4 5];
C.Zser  = [0.02+0.06i; 0.08+0.24i; 0.06+0.25i; 0.06+0.18i; 0.04+0.12i; 0.01+0.03i; 0.08+0.24i];
C.Bend  = [0.03i;      0.025i;     0.020i;     0.020i;     0.015i;     0.010i;     0.025i];

% Base loads (MW / MVAr) – generators other than PV/wind are 0,
% slack will balance. (PV on bus 5, wind on bus 4 added later)
C.Pd_MW   = [0 20 45 40 60];
C.Qd_MVAr = [0 10 15  5 10];
C.Pg_MW   = [0  0  0  0  0];
C.Qg_MVAr = [0 30  0  0  0];   % only bus-2 has fixed +30 MVAr

C.svc_buses = 2:5;
C.mid_list  = C.pairs;


here = fileparts(mfilename('fullpath'));
projectRoot = fileparts(here);                 % parent of code/
C.data.pv   = fullfile(projectRoot,'data','pvdata2.csv');
C.data.wind = fullfile(projectRoot,'data','windspeed.csv');

%  PV module (datasheet) 
C.PV.N_modules = round(30e6 / 280);           % 30 MW farm, 280 W per module
C.PV.VMPP  = 31.28;
C.PV.IMPP  = 8.95;
C.PV.VOC   = 37.82;
C.PV.ISC   = 9.47;
C.PV.NOT   = 43;                               % °C
C.PV.Kv    = -0.0032 * C.PV.VOC;               % V/°C  (-0.32% of Voc per °C)
C.PV.Ki    = -0.0006 * C.PV.ISC;               % A/°C  (-0.06% of Isc per °C)

%  Wind farm (power curve; 15 MW rated)
C.WF.PrMW  = 15;
C.WF.v_in  = 3;           % m/s (cut-in)
C.WF.v_r   = 12.5;        % m/s (rated ~12–13)
C.WF.v_out = 25;          % m/s (cut-out)
end
