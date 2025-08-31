% NR Load Flow with an SVC at the midpoint of line 3-4 (new bus 6)
% Base: Sbase = 100 MVA

clc; clear; close all;


y12 = 1/(0.02 + 0.06i);   B12 = 0.03i;
y13 = 1/(0.08 + 0.24i);   B13 = 0.025i;
y23 = 1/(0.06 + 0.25i);   B23 = 0.020i;
y24 = 1/(0.06 + 0.18i);   B24 = 0.020i;
y25 = 1/(0.04 + 0.12i);   B25 = 0.015i;
y34 = 1/(0.01 + 0.03i);   B34 = 0.010i;   % per-end shunt susceptance at buses 3 & 4
y45 = 1/(0.08 + 0.24i);   B45 = 0.025i;

% Original 5-bus Ybus
Y5 = [ (y12+y13+B12+B13)                     -y12                     -y13                   0                      0
       -y12               (y12+y23+y24+y25+B12+B23+B24+B25)           -y23                 -y24                  -y25
       -y13                                 -y23        (y13+y23+y34+B13+B23+B34)         -y34                   0
        0                                   -y24                               -y34   (y34+y45+y24+B34+B45+B24) -y45
        0                                   -y25                                0                   -y45   (y25+y45+B25+B45) ];


Y = zeros(6,6);      % 6-bus network after adding the midpoint
Y(1:5,1:5) = Y5;

% remove original 3-4 series and its two end-shunts (they will be re-added as half-lines)
Y(3,3) = Y(3,3) - y34 - B34;
Y(4,4) = Y(4,4) - y34 - B34;
Y(3,4) = 0;  Y(4,3) = 0;

% two half-π sections: 3–6 and 6–4
z34    = 0.01 + 0.03i;
z_half = z34/2;
y_half = 1/z_half;        % = 2*y34
B_half = B34/2;           % per-end shunt on each half-line

% connect 3–6
Y(3,3) = Y(3,3) + y_half + B_half;
Y(6,6) = Y(6,6) + y_half + B_half;
Y(3,6) = Y(3,6) - y_half;   Y(6,3) = Y(6,3) - y_half;

% connect 6–4
Y(4,4) = Y(4,4) + y_half + B_half;
Y(6,6) = Y(6,6) + y_half + B_half;
Y(4,6) = Y(4,6) - y_half;   Y(6,4) = Y(6,4) - y_half;

N = 6;


Sbase   = 100;           % MVA
svc_bus = 6;             % midpoint bus (new)
Vref    = 1.00;          % |V| target at SVC bus
B_svc   = 0.00;          % initial susceptance
Bmin    = -1.00;         % ~ -100 MVAr
Bmax    = +1.00;         % ~ +100 MVAr

% Scheduled load / gen (bus 2 is PQ gen)
Pd = [0 20 45 40 60 0] / Sbase;   % pu
Qd = [0 10 15  5 10 0] / Sbase;
Pg = [0 40  0  0  0 0] / Sbase;
Qg = [0 30  0  0  0 0] / Sbase;

Psch = Pg - Pd;
Qsch = Qg - Qd;

% Flat start
Vmag = ones(1,N); Vmag(1)=1.06; Vmag(2:5)=1.00; Vmag(6)=Vref;
Vang = zeros(1,N);


tol = 1e-6; max_iter = 50; iter = 0; err = 1;
ang_idx = 2:N; PQ = N-1;

while err > tol && iter < max_iter
    iter = iter + 1;

    % Y including SVC susceptance at its bus (shunt to ground)
    Yit = Y;
    Yit(svc_bus,svc_bus) = Yit(svc_bus,svc_bus) + 1i*B_svc;

    % enforce |V6|
    Vmag(svc_bus) = Vref;

    % injections from current state
    P = zeros(1,N); Q = zeros(1,N);
    for i=1:N
        for k=1:N
            th = angle(Yit(i,k)) + Vang(k) - Vang(i);
            P(i) = P(i) + Vmag(i)*Vmag(k)*abs(Yit(i,k))*cos(th);
            Q(i) = Q(i) - Vmag(i)*Vmag(k)*abs(Yit(i,k))*sin(th);
        end
    end

    % mismatch  (buses 2..N)
    dP = Psch - P;  dQ = Qsch - Q;
    M  = [dP(ang_idx) dQ(ang_idx)]';
    err = max(abs(M)); if err <= tol, break; end

    % Jacobian blocks
    J1 = zeros(PQ,PQ); J2 = zeros(PQ,PQ);
    J3 = zeros(PQ,PQ); J4 = zeros(PQ,PQ);
    for i=2:N
        for k=2:N
            if i==k
                J1(i-1,k-1) = -Q(i) - Vmag(i)^2*abs(Yit(i,i))*sin(angle(Yit(i,i)));
                J2(i-1,k-1) =  P(i) + Vmag(i)^2*abs(Yit(i,i))*cos(angle(Yit(i,i)));
                J3(i-1,k-1) =  P(i) - Vmag(i)^2*abs(Yit(i,i))*cos(angle(Yit(i,i)));
                J4(i-1,k-1) =  Q(i) - Vmag(i)^2*abs(Yit(i,i))*sin(angle(Yit(i,i)));
            else
                th = angle(Yit(i,k)) + Vang(k) - Vang(i);
                J1(i-1,k-1) = -Vmag(i)*Vmag(k)*abs(Yit(i,k))*sin(th);
                J2(i-1,k-1) =  Vmag(i)*Vmag(k)*abs(Yit(i,k))*cos(th);
                J3(i-1,k-1) = -Vmag(i)*Vmag(k)*abs(Yit(i,k))*cos(th);
                J4(i-1,k-1) = -Vmag(i)*Vmag(k)*abs(Yit(i,k))*sin(th);
            end
        end
    end

    % Replace |V6| column by dB column (dP/dB=0, dQ/dB=-|V6|^2 on the SVC row)
    J2(:, svc_bus-1) = 0;
    J4(:, svc_bus-1) = 0;
    J4(svc_bus-1, svc_bus-1) = -Vmag(svc_bus)^2;

    % solve & update
    dx = ( [J1 J2; J3 J4] \ M ).';
    Vang(ang_idx) = Vang(ang_idx) + dx(1:PQ);

    dVm = dx(PQ+1:end);
    if svc_bus > 2, Vmag(2:svc_bus-1) = Vmag(2:svc_bus-1) + dVm(1:svc_bus-2); end
    if svc_bus < N, Vmag(svc_bus+1:N) = Vmag(svc_bus+1:N) + dVm(svc_bus:end); end

    dB    = dVm(svc_bus-1);
    B_svc = max(min(B_svc + dB, Bmax), Bmin);
end


Yfinal = Y; Yfinal(svc_bus,svc_bus) = Yfinal(svc_bus,svc_bus) + 1i*B_svc;
Vdeg   = Vang*180/pi;

fprintf('Converged in %d iterations. Max mismatch = %.3e\n', iter, err);
for b=1:N
    fprintf('Bus %d: |V| = %.5f  angle = %+8.4f deg\n', b, Vmag(b), Vdeg(b));
end
Qsvc = -(Vmag(svc_bus)^2)*B_svc*Sbase;    % MVAr, + = capacitive injection
fprintf('\nSVC @ bus %d: B = %+8.5f pu   =>   Q_svc = %+8.3f MVAr\n\n', svc_bus, B_svc, Qsvc);


lines = [1 2; 1 3; 2 3; 2 4; 2 5; 3 6; 6 4; 4 5];
Zs    = [0.02+0.06i; 0.08+0.24i; 0.06+0.25i; 0.06+0.18i; ...
         0.04+0.12i; z_half;     z_half;     0.08+0.24i];

L = size(lines,1);
Iabs = zeros(L,1); Iang = zeros(L,1);
Pl   = zeros(L,1); Ql   = zeros(L,1);
for e=1:L
    i = lines(e,1); k = lines(e,2);
    [Iabs(e),Iang(e),Pl(e),Ql(e)] = current_and_lineloss(Vmag(i),Vdeg(i),Vmag(k),Vdeg(k),Zs(e),Sbase);
end
totPloss = sum(Pl);  totQloss = sum(Ql);


Pg1 = Sbase * calc_injection_P(1, Vmag, Vang, Yfinal);
Qg1 = Sbase * calc_injection_Q(1, Vmag, Vang, Yfinal);


% Table A: Bus results (generation & voltage)
genMW   = zeros(N,1); genMVAr = zeros(N,1);
genMW(1)= Pg1; genMVAr(1)= Qg1;
genMW(2)= 40;  genMVAr(2)= 30;   % scheduled PQ gen at bus 2

busNo   = (1:N).';
busVolt = Vmag(:);
busAng  = Vdeg(:);

T_A = table(busNo, genMW, genMVAr, busVolt, busAng, ...
    'VariableNames', {'Bus','GenMW','GenMVAr','Voltage_pu','Angle_deg'});

fprintf('\n================ Table A: Load Flow with SVC (Variable Susceptance Model) ================\n');
disp(T_A)

% Table B: Branch currents & losses
frm = lines(:,1); to = lines(:,2);
T_B = table(frm, to, Iabs, Iang, Pl, Ql, ...
    'VariableNames', {'From','To','Current_pu','CurrentAngle_deg','RealLoss_MW','ReactiveLoss_MVAr'});

fprintf('\n================ Table B: Current in each line and line loss ================\n');
disp(T_B)

% Table C: Per-bus totals
% (Allocate each line's loss half to each incident bus; this is a common, symmetric convention.)
lossMW_bus  = zeros(N,1);
lossMVAr_bus= zeros(N,1);
for e=1:L
    i = lines(e,1); k = lines(e,2);
    lossMW_bus([i k])   = lossMW_bus([i k])   + Pl(e)/2;
    lossMVAr_bus([i k]) = lossMVAr_bus([i k]) + Ql(e)/2;
end

loadMW    = Pd(:)*Sbase;
loadMVAr  = Qd(:)*Sbase;

T_C = table(busNo, genMW, genMVAr, loadMW, loadMVAr, lossMW_bus, lossMVAr_bus, ...
    'VariableNames', {'Bus','TotalGen_MW','TotalGen_MVAr','TotalLoad_MW','TotalLoad_MVAr','TotalLoss_MW','TotalLoss_MVAr'});

fprintf('\n================ Table C: Total load, total loss and total generation (per bus) ================\n');
disp(T_C)

% System totals & SVC rating (pretty print)
sysGenMW   = sum(genMW);
sysGenMVAr = sum(genMVAr);
sysLoadMW  = sum(loadMW);
sysLoadMVAr= sum(loadMVAr);
sysLossMW  = totPloss;     % same as sum(lossMW_bus)
sysLossMVAr= totQloss;

fprintf('---------------- System totals ----------------\n');
fprintf('Generation = %.4f MW / %.4f MVAr\n', sysGenMW, sysGenMVAr);
fprintf('Load       = %.4f MW / %.4f MVAr\n', sysLoadMW, sysLoadMVAr);
fprintf('Line loss  = %.4f MW / %.4f MVAr\n', sysLossMW, sysLossMVAr);

svcType = 'Inductive'; if Qsvc >= 0, svcType = 'Capacitive'; end
fprintf('\nSVC rating = %.4f MVAr (%s)\n\n', abs(Qsvc), svcType);


function [Iabs,Iang,Ploss,Qloss] = current_and_lineloss(Vm_i,ang_i_deg,Vm_j,ang_j_deg,Z,Sbase)
    Vi = Vm_i * exp(1i*deg2rad(ang_i_deg));
    Vj = Vm_j * exp(1i*deg2rad(ang_j_deg));
    I  = (Vi - Vj) / Z;
    Iabs = abs(I); Iang = rad2deg(angle(I));
    S_loss_pu = (Iabs^2) * Z;     % I^2 * Z (series loss)
    Ploss = Sbase * real(S_loss_pu);
    Qloss = Sbase * imag(S_loss_pu);
end

function P = calc_injection_P(i, Vmag, Vang, Y)
    N = numel(Vmag); P = 0;
    for k = 1:N
        th = angle(Y(i,k)) + Vang(k) - Vang(i);
        P = P + Vmag(i)*Vmag(k)*abs(Y(i,k))*cos(th);
    end
end

function Q = calc_injection_Q(i, Vmag, Vang, Y)
    N = numel(Vmag); Q = 0;
    for k = 1:N
        th = angle(Y(i,k)) + Vang(k) - Vang(i);
        Q  = Q - Vmag(i)*Vmag(k)*abs(Y(i,k))*sin(th);
    end
end
