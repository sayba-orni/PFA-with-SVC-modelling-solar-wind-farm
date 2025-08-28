% NR Load Flow with SVC at Bus 3 (Bus 2 is PQ) 
% Sbase = 100 MVA

clc; clear; close all;

y12 = 1/(0.02 + 0.06i);   B12 = 0.03i;
y13 = 1/(0.08 + 0.24i);   B13 = 0.025i;
y23 = 1/(0.06 + 0.25i);   B23 = 0.020i;
y24 = 1/(0.06 + 0.18i);   B24 = 0.020i;
y25 = 1/(0.04 + 0.12i);   B25 = 0.015i;
y34 = 1/(0.01 + 0.03i);   B34 = 0.010i;
y45 = 1/(0.08 + 0.24i);   B45 = 0.025i;

Y_base = [ (y12+y13+B12+B13)                     -y12                     -y13                   0                      0
           -y12               (y12+y23+y24+y25+B12+B23+B24+B25)           -y23                 -y24                  -y25
           -y13                                 -y23        (y13+y23+y34+B13+B23+B34)         -y34                   0
            0                                   -y24                               -y34   (y34+y45+y24+B34+B45+B24) -y45
            0                                   -y25                                0                   -y45   (y25+y45+B25+B45) ];

Sbase     = 100;                % MVA
slack     = 1;
svc_bus   = 3;                  % SVC shunt at bus 3
Vref_svc  = 1.00;               % target |V3| while regulating
B_svc     = 0.0;                % initial susceptance
Bmin      = -1.0;               % ≈ -100 MVAr   (capacitive => injects +Q)
Bmax      = +1.0;               % ≈ +100 MVAr   (inductive  => absorbs Q)
svc_regulating = true;          % start regulating

% Flat start
V_abs   = [1.06, 1.00, Vref_svc, 1.00, 1.00];
V_ang   = zeros(1,5);

% Scheduled load/gen (bus 2 is PQ here)
Pd = [0  20 45 40 60]/Sbase;
Qd = [0  10 15  5 10]/Sbase;
Pg = [0  40  0  0  0]/Sbase;
Qg = [0  30  0  0  0]/Sbase;     % Q at bus 2 is fixed (PQ) in this variant

Psch = Pg - Pd;
Qsch = Qg - Qd;

%%  NR settings 
tol = 1e-6; max_iter = 50; iter = 0; err = 1;

ang_idx = 2:5;                 % angles for buses 2..5
pq_set  = 2:5;                 % PQ buses (all non-slack)

while err > tol && iter < max_iter
    iter = iter + 1;

    %  Ybus with SVC this iteration 
    Y = Y_base;
    Y(svc_bus, svc_bus) = Y(svc_bus, svc_bus) + 1i*B_svc;

    %  Compute P, Q from current state 
    N = 5;
    Pcalc = zeros(1,N); Qcalc = zeros(1,N);
    for i = 1:N
        for k = 1:N
            Vi = V_abs(i); Vk = V_abs(k);
            Yik = Y(i,k); th = angle(Yik) + V_ang(k) - V_ang(i);
            Pcalc(i) = Pcalc(i) + Vi*Vk*abs(Yik)*cos(th);
            Qcalc(i) = Qcalc(i) - Vi*Vk*abs(Yik)*sin(th);
        end
    end

    % Hold |V3| to setpoint while regulating
    if svc_regulating, V_abs(svc_bus) = Vref_svc; end

    %  Mismatch vector 
    dP = Psch - Pcalc;
    dQ = Qsch - Qcalc;
    misP = dP(ang_idx);          % ΔP for 2..5
    misQ = dQ(pq_set);           % ΔQ for 2..5
    M = [misP, misQ]';           % 8x1

    err = max(abs(M)); if err <= tol, break; end

    % Build Jacobian 
    na = numel(ang_idx);               % 4
    mag_unknowns = pq_set;             % start with 2..5 magnitudes
    addB = 0;
    if svc_regulating
        % replace |V3| column by dB column
        mag_unknowns = setdiff(mag_unknowns, svc_bus);
        addB = 1;
    end
    nm = numel(mag_unknowns);

    J11 = zeros(na,na);                % dP/dδ
    J12 = zeros(na,nm+addB);           % dP/d|V| (+ dP/dB column if addB)
    J21 = zeros(numel(pq_set), na);    % dQ/dδ
    J22 = zeros(numel(pq_set), nm+addB);% dQ/d|V| (+ dQ/dB column if addB)

    % Fill J11, J12
    for r = 1:na
        i = ang_idx(r);
        for c = 1:na
            k = ang_idx(c);
            if i==k
                J11(r,c) = -Qcalc(i) - V_abs(i)^2*abs(Y(i,i))*sin(angle(Y(i,i)));
            else
                thik = angle(Y(i,k)) + V_ang(k) - V_ang(i);
                J11(r,c) = -V_abs(i)*V_abs(k)*abs(Y(i,k))*sin(thik);
            end
        end
        for c = 1:nm
            k = mag_unknowns(c);
            if i==k
                J12(r,c) =  Pcalc(i) + V_abs(i)^2*abs(Y(i,i))*cos(angle(Y(i,i)));
            else
                thik = angle(Y(i,k)) + V_ang(k) - V_ang(i);
                J12(r,c) =  V_abs(i)*V_abs(k)*abs(Y(i,k))*cos(thik);
            end
        end
    end
    % dP/dB column is zero
    if addB, J12(:, nm+1) = 0; end

    % Fill J21, J22 (rows correspond to pq_set)
    for r = 1:numel(pq_set)
        i = pq_set(r);
        for c = 1:na
            k = ang_idx(c);
            if i==k
                J21(r,c) =  Pcalc(i) - V_abs(i)^2*abs(Y(i,i))*cos(angle(Y(i,i)));
            else
                thik = angle(Y(i,k)) + V_ang(k) - V_ang(i);
                J21(r,c) = -V_abs(i)*V_abs(k)*abs(Y(i,k))*cos(thik);
            end
        end
        for c = 1:nm
            k = mag_unknowns(c);
            if i==k
                J22(r,c) =  Qcalc(i) - V_abs(i)^2*abs(Y(i,i))*sin(angle(Y(i,i)));
            else
                thik = angle(Y(i,k)) + V_ang(k) - V_ang(i);
                J22(r,c) = -V_abs(i)*V_abs(k)*abs(Y(i,k))*sin(thik);
            end
        end
    end

    % dQ/dB column: only the SVC bus row is +|V3|^2
    if addB
        J22(:, nm+1) = 0;
        row_svc = find(pq_set==svc_bus);
        J22(row_svc, nm+1) = V_abs(svc_bus)^2;
    end

    % Solve
    J = [J11 J12; J21 J22];
    dx = J \ M;

    % Updates
    dDel    = dx(1:na).';
    dRest   = dx(na+1:end).';

    V_ang(ang_idx) = V_ang(ang_idx) + dDel;

    if nm > 0
        dVrel = dRest(1:nm);
        V_abs(mag_unknowns) = V_abs(mag_unknowns) .* (1 + dVrel);
    end
    if addB
        dB = dRest(end);
        B_svc = B_svc + dB;
        % If B hits limit, clamp and release regulation next iteration
        if B_svc > Bmax, B_svc = Bmax; svc_regulating = false; end
        if B_svc < Bmin, B_svc = Bmin; svc_regulating = false; end
    end

    % enforce the target while regulating
    if svc_regulating, V_abs(svc_bus) = Vref_svc; end
end

%%  Results
V_deg = V_ang * 180/pi;
fprintf('Converged in %d iterations. Max mismatch = %.3e\n', iter, err);
for b = 1:5
    fprintf('Bus %d: |V| = %.5f  angle = %+8.4f deg\n', b, V_abs(b), V_deg(b));
end
Qsvc_MVAr = -(V_abs(svc_bus)^2) * B_svc * Sbase;
fprintf('\nSVC: B = %+8.5f pu   =>   Q_svc = %+8.3f MVAr\n', B_svc, Qsvc_MVAr);

%%  line currents & losses (series only) 
lines = [1 2; 1 3; 2 3; 2 4; 2 5; 3 4; 4 5];
Zs    = [0.02+0.06i; 0.08+0.24i; 0.06+0.25i; 0.06+0.18i; ...
         0.04+0.12i; 0.01+0.03i; 0.08+0.24i];

fprintf('\nFrom To\t I_mag(pu)\t I_ang(deg)\t P_loss(MW)\t Q_loss(MVAr)\n');
for e = 1:size(lines,1)
    i = lines(e,1); k = lines(e,2);
    out = current_and_lineloss(V_abs(i), V_deg(i), V_abs(k), V_deg(k), Zs(e), Sbase);
    fprintf('%d    %d\t %.5f\t  %8.4f\t   %9.5f\t   %9.5f\n', i, k, out(1), out(2), out(3), out(4));
end

%  helper 
function r = current_and_lineloss(Vm_i, ang_i_deg, Vm_j, ang_j_deg, Z, Sbase)
    Vi = Vm_i * exp(1i*deg2rad(ang_i_deg));
    Vj = Vm_j * exp(1i*deg2rad(ang_j_deg));
    Iij = (Vi - Vj) / Z;                 % current i->j (series branch)
    Iabs = abs(Iij); Iang = rad2deg(angle(Iij));
    S_loss_pu = (Iabs^2) * Z;            % I^2*Z
    Ploss_MW  = Sbase * real(S_loss_pu);
    Qloss_MVAr= Sbase * imag(S_loss_pu);
    r = [Iabs, Iang, Ploss_MW, Qloss_MVAr];
end

