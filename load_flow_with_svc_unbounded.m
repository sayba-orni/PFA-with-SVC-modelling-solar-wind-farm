clc;
clear all;
close all;


y12=(0.02+0.06*i)^-1;  B12=0.03*i;
y13=(0.08+0.24*i)^-1;  B13=0.025*i;
y23=(0.06+0.25*i)^-1;  B23=0.02*i;
y24=(0.06+0.18*i)^-1;  B24=0.02*i;
y25=(0.04+0.12*i)^-1;  B25=0.015*i;
y34=(0.01+0.03*i)^-1;  B34=0.01*i;
y45=(0.08+0.24*i)^-1;  B45=0.025*i;

Y_base=[ y12+y13+B12+B13   -y12                       -y13                 0                       0 ;
        -y12         y12+y23+y24+y25+B23+B25+B24+B12  -y23               -y24                    -y25
        -y13               -y23            y13+y23+y34+B13+B23+B34       -y34                     0
         0                 -y24                      -y34          y34+y45+y24+B24+B34+B45       -y45
         0                 -y25                       0                  -y45              y25+y45+B25+B45];


svc_bus  = 3;
Vref_svc = 1.00;         % enforce |V3|
B_svc    = 1;            % initial susceptance (pu)


V_abs   = [1.06 1 1 1 1];
V_angle = [0 0 0 0 0];
Sbase   = 100;           % MVA


Pd_pu = [0 20 45 40 60]/Sbase;
Qd_pu = [0 10 15  5 10]/Sbase;

Pg_pu = [0 40  0  0  0]/Sbase;   % slack Pg is unknown; others fixed
Qg_pu = [0 30  0  0  0]/Sbase;   % bus 2 is PQ in this variant

P_scheduled = Pg_pu - Pd_pu;
Q_scheduled = Qg_pu - Qd_pu;


tol=1e-6; num_of_it=1; err=1;

n=5; PV=0; PQ=n-PV-1;
J1=zeros(PV+PQ,PV+PQ);
J2=zeros(PV+PQ,PQ);
J3=zeros(PQ,PV+PQ);
J4=zeros(PQ,PQ);

while err>=tol
    % Ybus with the current SVC B
    Y_matrix = Y_base;
    Y_matrix(svc_bus, svc_bus) = Y_base(svc_bus, svc_bus) + 1i*B_svc;

    % hold |V3| at Vref during regulation
    V_abs(svc_bus) = Vref_svc;

    % P,Q injections from present state
    P_calculated=zeros(1,5); Q_calculated=zeros(1,5);
    for ii=1:5
        for kk=1:5
            th = angle(Y_matrix(ii,kk)) + V_angle(kk) - V_angle(ii);
            P_calculated(ii) = P_calculated(ii) + V_abs(ii)*V_abs(kk)*abs(Y_matrix(ii,kk))*cos(th);
            Q_calculated(ii) = Q_calculated(ii) - V_abs(ii)*V_abs(kk)*abs(Y_matrix(ii,kk))*sin(th);
        end
    end

    % mismatch (buses 2..5 only)
    del_P = P_scheduled - P_calculated;
    del_Q = Q_scheduled - Q_calculated;
    Mismatch = [del_P(2:5), del_Q(2:5)]';
    err = max(abs(Mismatch));
    if err<tol, break; end

    % Jacobians
    % J1 = dP/dδ, J2 = dP/d|V|, J3 = dQ/dδ, J4 = dQ/d|V|
    for ii = 2:5
        for jj = 2:5
            if ii==jj
                J1(ii-1,jj-1) = -Q_calculated(ii) - V_abs(ii)^2*abs(Y_matrix(ii,ii))*sin(angle(Y_matrix(ii,ii)));
                J2(ii-1,jj-1) =  P_calculated(ii) + V_abs(ii)^2*abs(Y_matrix(ii,ii))*cos(angle(Y_matrix(ii,ii)));
                J3(ii-1,jj-1) =  P_calculated(ii) - V_abs(ii)^2*abs(Y_matrix(ii,ii))*cos(angle(Y_matrix(ii,ii)));
                J4(ii-1,jj-1) =  Q_calculated(ii) - V_abs(ii)^2*abs(Y_matrix(ii,ii))*sin(angle(Y_matrix(ii,ii)));
            else
                th = angle(Y_matrix(ii,jj)) + V_angle(jj) - V_angle(ii);
                J1(ii-1,jj-1) = -V_abs(ii)*V_abs(jj)*abs(Y_matrix(ii,jj))*sin(th);
                J2(ii-1,jj-1) =  V_abs(ii)*V_abs(jj)*abs(Y_matrix(ii,jj))*cos(th);
                J3(ii-1,jj-1) = -V_abs(ii)*V_abs(jj)*abs(Y_matrix(ii,jj))*cos(th);
                J4(ii-1,jj-1) = -V_abs(ii)*V_abs(jj)*abs(Y_matrix(ii,jj))*sin(th);
            end
        end
    end
    
    J2(:, svc_bus-1) = 0;                         % dP/dB = 0
    J4(:, svc_bus-1) = 0;  J4(svc_bus-1,svc_bus-1) = -V_abs(svc_bus)^2;

   
    J = [J1 J2; J3 J4];
    dx = J \ Mismatch;
    dx = dx.';                  % row vector

    % angles
    V_angle(2:5) = V_angle(2:5) + dx(1:4);

    % magnitudes 
    V_abs(2:svc_bus-1)   = V_abs(2:svc_bus-1)   + dx(5:(4+svc_bus-2));
    V_abs(svc_bus+1:end) = V_abs(svc_bus+1:end) + dx((4+svc_bus):end);

    % susceptance update
    dB   = dx(4 + svc_bus - 1);
    B_svc = B_svc + dB;

    num_of_it = num_of_it + 1;
end



Y_matrix = Y_base;  Y_matrix(svc_bus,svc_bus)=Y_base(svc_bus,svc_bus)+1i*B_svc;
P_calculated=zeros(1,5); Q_calculated=zeros(1,5);
for ii=1:5
    for kk=1:5
        th = angle(Y_matrix(ii,kk)) + V_angle(kk) - V_angle(ii);
        P_calculated(ii) = P_calculated(ii) + V_abs(ii)*V_abs(kk)*abs(Y_matrix(ii,kk))*cos(th);
        Q_calculated(ii) = Q_calculated(ii) - V_abs(ii)*V_abs(kk)*abs(Y_matrix(ii,kk))*sin(th);
    end
end
iter = num_of_it-1;

% bus volt/angle
V_deg = V_angle*180/pi;

Pg1_MW  =  Sbase * P_calculated(1);
Qg1_MVAr=  Sbase * Q_calculated(1);

Qsvc_MVAr = -(V_abs(svc_bus)^2) * B_svc * Sbase;

fprintf('Converged in %d iterations. Max mismatch = %.3e\n', iter, err);
for b=1:5
    fprintf('Bus %d: |V| = %.5f  angle = %8.4f deg\n', b, V_abs(b), V_deg(b));
end
fprintf('\nSVC @ bus %d: B = %+8.5f pu   =>   Q_svc = %+8.3f MVAr\n\n', svc_bus, B_svc, Qsvc_MVAr);

Current_and_lineloss=zeros(7,4);
Current_and_lineloss(1,:)=current_and_lineloss(V_abs(1),V_deg(1),V_abs(2),V_deg(2),0.02,0.06);
Current_and_lineloss(2,:)=current_and_lineloss(V_abs(1),V_deg(1),V_abs(3),V_deg(3),0.08,0.24);
Current_and_lineloss(3,:)=current_and_lineloss(V_abs(2),V_deg(2),V_abs(3),V_deg(3),0.06,0.25);
Current_and_lineloss(4,:)=current_and_lineloss(V_abs(2),V_deg(2),V_abs(4),V_deg(4),0.06,0.18);
Current_and_lineloss(5,:)=current_and_lineloss(V_abs(2),V_deg(2),V_abs(5),V_deg(5),0.04,0.12);
Current_and_lineloss(6,:)=current_and_lineloss(V_abs(3),V_deg(3),V_abs(4),V_deg(4),0.01,0.03);
Current_and_lineloss(7,:)=current_and_lineloss(V_abs(4),V_deg(4),V_abs(5),V_deg(5),0.08,0.24);

fprintf('From To\t I_mag(pu)\t I_ang(deg)\t P_loss(MW)\t Q_loss(MVAr)\n');
lines={'1  2','1  3','2  3','2  4','2  5','3  4','4  5'};
for k=1:size(Current_and_lineloss,1)
    r=Current_and_lineloss(k,:);
    fprintf('%s\t %.5f\t %10.4f\t %10.5f\t %10.5f\n',lines{k},r(1),r(2),r(3),r(4));
end
totPloss = sum(Current_and_lineloss(:,3));
totQloss = sum(Current_and_lineloss(:,4));
fprintf('\nTotal line loss:  P = %.5f MW,   Q = %.5f MVAr\n\n', totPloss, totQloss);


% Generation (MW/MVAr)
genMW   = [Pg1_MW, 40, 0, 0, 0];
genMVAr = [Qg1_MVAr, 30, 0, 0, 0];

% Load (MW/MVAr)
loadMW   = [0, 20, 45, 40, 60];
loadMVAr = [0, 10, 15,  5, 10];

fprintf('Bus   Generation(MW,MVAr)     Voltage(pu)   Angle(deg)\n');
for b=1:5
    fprintf('%-4d  %8.4f  %8.4f      %8.4f      %8.4f\n', ...
        b, genMW(b), genMVAr(b), V_abs(b), V_deg(b));
end
fprintf('\nTotals: Gen = %.4f MW / %.4f MVAr,  Load = %.4f MW / %.4f MVAr,  Line Loss = %.4f MW / %.4f MVAr\n', ...
    sum(genMW), sum(genMVAr), sum(loadMW), sum(loadMVAr), totPloss, totQloss);

Lmw  = Current_and_lineloss(:,3);
Lmvar= Current_and_lineloss(:,4);


bus_lbl   = {'1(swing)','2','3','4','5','',''};
genMW_row = [Pg1_MW, 40, 0, 0, 0, 0, 0];
genMV_row = [Qg1_MVAr, 30, 0, 0, 0, 0, 0];

loadMW_row  = [0, 20, 45, 40, 60, 0, 0];
loadMV_row  = [0, 10, 15,  5, 10, 0, 0];

lossMW_row  = [Lmw(1), Lmw(2), Lmw(3), Lmw(4), Lmw(5), Lmw(6), Lmw(7)];
lossMV_row  = [Lmvar(1), Lmvar(2), Lmvar(3), Lmvar(4), Lmvar(5), Lmvar(6), Lmvar(7)];

fprintf('\nTotal load, total loss and total generation is shown in the following table:\n\n');
fprintf('%-10s %-20s %-20s %-20s\n','Bus No','Total Generation','Total Load','Total Loss(Line)');
fprintf('%-10s %-9s %-10s %-9s %-10s %-9s %-10s\n','','MW','MVAR','MW','MVAR','MW','MVAR');

for r = 1:numel(bus_lbl)
    fprintf('%-10s %9.4f %10.4f %9.0f %10.0f %9.4f %10.4f\n', ...
        bus_lbl{r}, genMW_row(r), genMV_row(r), ...
        loadMW_row(r), loadMV_row(r), lossMW_row(r), lossMV_row(r));
end


totGenMW   = sum(genMW_row);
totGenMVAr = sum(genMV_row);
totLoadMW  = sum(loadMW_row);
totLoadMVAr= sum(loadMV_row);
totLossMW  = sum(lossMW_row);
totLossMVAr= sum(lossMV_row);

fprintf('\nTotals: Gen = %.4f MW / %.4f MVAr, Load = %.4f MW / %.4f MVAr, Line Loss = %.4f MW / %.4f MVAr\n', ...
    totGenMW, totGenMVAr, totLoadMW, totLoadMVAr, totLossMW, totLossMVAr);


svc_abs   = abs(Qsvc_MVAr);
svc_type  = ternary(Qsvc_MVAr >= 0, 'Inductive', 'Capacitive'); % +Q = absorbs (inductive)
fprintf('\nSvc Rating = %.4f MVAR (%s)\n', svc_abs, svc_type);


function out = ternary(cond, a, b), if cond, out=a; else, out=b; end, end
