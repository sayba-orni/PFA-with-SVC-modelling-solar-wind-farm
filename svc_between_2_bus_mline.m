% NR Load Flow with an SVC at the midpoint of line 3-4 (new bus 6)
% Base: Sbase = 100 MVA

clc; clear; close all;


y12 = 1/(0.02 + 0.06i);   B12 = 0.03i;
y13 = 1/(0.08 + 0.24i);   B13 = 0.025i;

y23 = 1/(0.06 + 0.25i);   B23 = 0.020i;
y24 = 1/(0.06 + 0.18i);   B24 = 0.020i;
y25 = 1/(0.04 + 0.12i);   B25 = 0.015i;
y34 = 1/(0.01 + 0.03i);   B34 = 0.010i;   % NOTE: B34 is the end-shunt used at bus 3 and 4
y45 = 1/(0.08 + 0.24i);   B45 = 0.025i;

% Original 5-bus Y
Y5 = [ (y12+y13+B12+B13)                     -y12                     -y13                   0                      0
       -y12               (y12+y23+y24+y25+B12+B23+B24+B25)           -y23                 -y24                  -y25
       -y13                                 -y23        (y13+y23+y34+B13+B23+B34)         -y34                   0
        0                                   -y24                               -y34   (y34+y45+y24+B34+B45+B24) -y45
        0                                   -y25                                0                   -y45   (y25+y45+B25+B45) ];


Y = zeros(6,6);
Y(1:5,1:5) = Y5;

% Remove the original 3-4 branch COMPLETELY (series AND the per-end shunts)
Y(3,3) = Y(3,3) - y34 - B34;
Y(4,4) = Y(4,4) - y34 - B34;
Y(3,4) = 0;  Y(4,3) = 0;

% Add two half-π sections: 3–6 and 6–4
z34     = 0.01 + 0.03i;
z_half  = z34/2;
y_half  = 1/z_half;       % = 2*y34
B_half  = B34/2;          % each half-line end-shunt

% 3–6
Y(3,3) = Y(3,3) + y_half + B_half;
Y(6,6) = Y(6,6) + y_half + B_half;
Y(3,6) = Y(3,6) - y_half;   Y(6,3) = Y(6,3) - y_half;

% 6–4
Y(4,4) = Y(4,4) + y_half + B_half;
Y(6,6) = Y(6,6) + y_half + B_half;
Y(4,6) = Y(4,6) - y_half;   Y(6,4) = Y(6,4) - y_half;

N = 6;                      % total buses after adding the midpoint

Sbase     = 100;            % MVA
slack     = 1;
svc_bus   = 6;              % the midpoint bus
Vref_svc  = 1.00;
B_svc     = 0.00;           % initial SVC susceptance
Bmin      = -1.00;          % ~ -100 MVAr
Bmax      = +1.00;          % ~ +100 MVAr

% Scheduled load/gen (bus 2 is PQ; others as in classic 5-bus case)
Pd = [0 20 45 40 60]/Sbase;     Pd = [Pd 0];   % no load at bus 6
Qd = [0 10 15  5 10]/Sbase;     Qd = [Qd 0];
Pg = [0 40  0  0  0]/Sbase;     Pg = [Pg 0];
Qg = [0 30  0  0  0]/Sbase;     Qg = [Qg 0];

Psch = Pg - Pd;     Qsch = Qg - Qd;

% Flat start
Vmag = ones(1,N);  Vmag(1)=1.06; Vmag(2:5)=1.00; Vmag(6)=Vref_svc;
Vang = zeros(1,N);

tol = 1e-6;  max_iter = 50;  iter = 0;  err = 1;
ang_idx = 2:N;      % unknown angles
PQ     = N-1;       % number of non-slack buses (angles and |V| columns)

while err > tol && iter < max_iter
    iter = iter + 1;

    % Y with SVC this iteration
    Yit = Y;
    Yit(svc_bus,svc_bus) = Yit(svc_bus,svc_bus) + 1i*B_svc;

    % hold |V6| = Vref
    Vmag(svc_bus) = Vref_svc;

    % injections
    Pcalc = zeros(1,N); Qcalc = zeros(1,N);
    for i = 1:N
        for k = 1:N
            th = angle(Yit(i,k)) + Vang(k) - Vang(i);
            Pcalc(i) = Pcalc(i) + Vmag(i)*Vmag(k)*abs(Yit(i,k))*cos(th);
            Qcalc(i) = Qcalc(i) - Vmag(i)*Vmag(k)*abs(Yit(i,k))*sin(th);
        end
    end

    % mismatches for buses 2..N
    dP = Psch - Pcalc;  dQ = Qsch - Qcalc;
    M  = [dP(ang_idx) dQ(ang_idx)]';
    err = max(abs(M));  if err<=tol, break; end

    % Jacobians (PQxPQ blocks)
    J1 = zeros(PQ,PQ); J2 = zeros(PQ,PQ);
    J3 = zeros(PQ,PQ); J4 = zeros(PQ,PQ);
    for i = 2:N
        for k = 2:N
            if i==k
                J1(i-1,k-1) = -Qcalc(i) - Vmag(i)^2*abs(Yit(i,i))*sin(angle(Yit(i,i)));
                J2(i-1,k-1) =  Pcalc(i) + Vmag(i)^2*abs(Yit(i,i))*cos(angle(Yit(i,i)));
                J3(i-1,k-1) =  Pcalc(i) - Vmag(i)^2*abs(Yit(i,i))*cos(angle(Yit(i,i)));
                J4(i-1,k-1) =  Qcalc(i) - Vmag(i)^2*abs(Yit(i,i))*sin(angle(Yit(i,i)));
            else
                th = angle(Yit(i,k)) + Vang(k) - Vang(i);
                J1(i-1,k-1) = -Vmag(i)*Vmag(k)*abs(Yit(i,k))*sin(th);
                J2(i-1,k-1) =  Vmag(i)*Vmag(k)*abs(Yit(i,k))*cos(th);
                J3(i-1,k-1) = -Vmag(i)*Vmag(k)*abs(Yit(i,k))*cos(th);
                J4(i-1,k-1) = -Vmag(i)*Vmag(k)*abs(Yit(i,k))*sin(th);
            end
        end
    end

    % Replace |V_svc| column by dB column (dP/dB=0; dQ/dB = -|V|^2 at SVC row)
    J2(:, svc_bus-1) = 0;
    J4(:, svc_bus-1) = 0;
    J4(svc_bus-1, svc_bus-1) = -Vmag(svc_bus)^2;

    % Solve and update
    dx = ( [J1 J2; J3 J4] \ M ).';
    Vang(ang_idx) = Vang(ang_idx) + dx(1:PQ);

    dVm = dx(PQ+1:end);
    if svc_bus > 2
        Vmag(2:svc_bus-1) = Vmag(2:svc_bus-1) + dVm(1:svc_bus-2);
    end
    if svc_bus < N
        Vmag(svc_bus+1:N) = Vmag(svc_bus+1:N) + dVm(svc_bus:end);
    end
    dB    = dVm(svc_bus-1);
    B_svc = max(min(B_svc + dB, Bmax), Bmin);
end


Yfinal = Y; Yfinal(svc_bus,svc_bus) = Yfinal(svc_bus,svc_bus) + 1i*B_svc;
Vdeg   = Vang*180/pi;

fprintf('Converged in %d iterations. Max mismatch = %.3e\n', iter, err);
for b=1:N
    fprintf('Bus %d: |V| = %.5f  angle = %+8.4f deg\n', b, Vmag(b), Vdeg(b));
end
Qsvc = -(Vmag(svc_bus)^2)*B_svc*Sbase;
fprintf('\nSVC @ bus %d: B = %+8.5f pu   =>   Q_svc = %+8.3f MVAr\n\n', svc_bus, B_svc, Qsvc);


lines = [1 2; 1 3; 2 3; 2 4; 2 5; 3 6; 6 4; 4 5];
Zs    = [0.02+0.06i; 0.08+0.24i; 0.06+0.25i; 0.06+0.18i; ...
         0.04+0.12i; z_half;     z_half;     0.08+0.24i];

fprintf('From To\t I_mag(pu)\t I_ang(deg)\t P_loss(MW)\t Q_loss(MVAr)\n');
lossPQ = zeros(size(lines,1),2);
for e=1:size(lines,1)
    i = lines(e,1); k = lines(e,2);
    [Iabs,Iang,Pl,Ql] = current_and_lineloss(Vmag(i),Vdeg(i),Vmag(k),Vdeg(k),Zs(e),Sbase);
    fprintf('%d    %d\t %.5f\t %10.4f\t %10.5f\t %10.5f\n', i,k,Iabs,Iang,Pl,Ql);
    lossPQ(e,:) = [Pl Ql];
end
fprintf('\nTotal line loss:  P = %.5f MW,   Q = %.5f MVAr\n\n', sum(lossPQ(:,1)), sum(lossPQ(:,2)));

% Slack injections (with final Y)
Pg1 = Sbase * calc_injection_P(1, Vmag, Vang, Yfinal);
Qg1 = Sbase * calc_injection_Q(1, Vmag, Vang, Yfinal);
fprintf('Slack injections:  P = %.4f MW,  Q = %.4f MVAr\n', Pg1, Qg1);


function [Iabs,Iang,Ploss,Qloss] = current_and_lineloss(Vm_i,ang_i_deg,Vm_j,ang_j_deg,Z,Sbase)
    Vi = Vm_i*exp(1i*deg2rad(ang_i_deg));
    Vj = Vm_j*exp(1i*deg2rad(ang_j_deg));
    I  = (Vi - Vj)/Z;
    Iabs = abs(I); Iang = rad2deg(angle(I));
    S_loss_pu = (Iabs^2)*Z;
    Ploss = Sbase*real(S_loss_pu);
    Qloss = Sbase*imag(S_loss_pu);
end

function P = calc_injection_P(i,Vmag,Vang,Y)
    N = numel(Vmag); P=0;
    for k=1:N
        th = angle(Y(i,k)) + Vang(k) - Vang(i);
        P = P + Vmag(i)*Vmag(k)*abs(Y(i,k))*cos(th);
    end
end

function Q = calc_injection_Q(i,Vmag,Vang,Y)
    N = numel(Vmag); Q=0;
    for k=1:N
        th = angle(Y(i,k)) + Vang(k) - Vang(i);
        Q = Q - Vmag(i)*Vmag(k)*abs(Y(i,k))*sin(th);
    end
end
