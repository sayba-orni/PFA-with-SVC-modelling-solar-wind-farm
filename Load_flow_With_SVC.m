clc;
clear all;
close all;

% Define Admittance Matrix
y12=(0.02+0.06*i)^-1;       B12=0.03*i;
y13=(0.08+0.24*i)^-1;       B13=0.025*i;
y23=(0.06+0.25*i)^-1;       B23=0.02*i;
y24=(0.06+0.18*i)^-1;       B24=0.02*i;
y25=(0.04+0.12*i)^-1;       B25=0.015*i;
y34=(0.01+0.03*i)^-1;       B34=0.01*i;
y45=(0.08+0.24*i)^-1;       B45=0.025*i;

Y_matrix=[ y12+y13+B12+B13   -y12                         -y13                0                        0 ;
    -y12              y12+y23+y24+y25+B23+B25+B24+B12   -y23              -y24                      -y25
    -y13                -y23                y13+y23+y34+B13+B23+B34   -y34                       0
    0                 -y24                        -y34             y34+y45+y24+B24+B34+B45    -y45
    0                 -y25                           0                -y45             y25+y45+B25+B45];

svc_bus=3;

% svc_suceptance
B_svc=1;              % 100 VAr Rating (initial)

% Define Voltage Magnitudes and Angles
V_abs = [1.06 , 1 , 1 , 1, 1 ] ;
V_angle = [0 0 0 0 0] ;


max_ittr=1;
num_of_it=1;
tol=0.000001;

n=5;
PV=0;
PQ=n-PV-1;

J1=zeros(PV+PQ,PV+PQ);
J2=zeros(PV+PQ,PQ);
J3=zeros(PQ,PV+PQ);
J4=zeros(PQ,PQ);
error=1;
%qerror=1;

while error>=tol

    Y_matrix(svc_bus,svc_bus)=Y_matrix(svc_bus,svc_bus)+B_svc*1i;

    % Define Power Demand 
    Pd_pu = [0  20  45  40 60 ].*1/100 ;
    Qd_pu = [0  10  15  5  10 ].*1/100 ;

    % Define Power Generation
    Pg_Pu = [0 , 40 , 0 , 0, 0]/100;

    Qg_pu = [0 , 0.30 , 0, 0, 0];
    if svc_bus==2

    Qg_pu(svc_bus)=0.3+(V_abs(svc_bus)*V_abs(svc_bus)*B_svc);
    else
      Qg_pu(svc_bus)=+V_abs(svc_bus)*V_abs(svc_bus)*B_svc;
    end

    P_scheduled = Pg_Pu - Pd_pu;
    Q_scheduled = Qg_pu - Qd_pu;


    % Calculate Calculated Power
    Q_calculated = zeros(1,5);
    P_calculated = zeros(1,5);

    for i = 1 : 5
        for k = 1 : 5
            Q_calculated(i) = Q_calculated(i) - V_abs(i)*V_abs(k)*abs(Y_matrix(i,k)) * sin( angle(Y_matrix(i,k)) + V_angle(k) - V_angle(i) );
            P_calculated(i) = P_calculated(i) + V_abs(i)*V_abs(k)*abs(Y_matrix(i,k)) * cos( angle(Y_matrix(i,k)) + V_angle(k) - V_angle(i) );
        end
    end
    % Calculate Mismatches
    del_P = P_scheduled - P_calculated;
    del_Q = Q_scheduled - Q_calculated;
    Mismatch = [del_P(2:5),del_Q(2:5)]';

    % Calculate Jacobians
    [r1,c1]=size(J1);
    for i = 2 : r1+1
        for j = 2 : c1+1
            if i == j
                J1(i-1,j-1) = - Q_calculated(i) - abs( V_abs(i)*V_abs(i)*abs(Y_matrix(i,i)) ) * sin( angle(Y_matrix(i,i)) );
            else
                J1(i-1,j-1) = -abs( V_abs(i)*V_abs(j)*abs(Y_matrix(i,j)) ) * sin( angle(Y_matrix(i,j)) + V_angle(j) - V_angle(i) );
            end
        end
    end

    [r1,c1]=size(J2);
    for i = 2 : r1+1
        for j = 2 : c1+1
            if i == j
                J2(i-1,j-1) = P_calculated(i) + abs( V_abs(i)*V_abs(i)*abs(Y_matrix(i,i)) ) * cos( angle(Y_matrix(i,i)) );
            else
                J2(i-1,j-1) = abs( V_abs(i)*V_abs(j)*abs(Y_matrix(i,j)) ) * cos( angle(Y_matrix(i,j)) + V_angle(j) - V_angle(i) );
            end
        end
    end
    J2_col=zeros(r1,1);
    J2(:,svc_bus-1)=J2_col;
    

    [r1,c1]=size(J3);
    for i = 2 : r1+1
        for j = 2 : c1+1
            if i == j
                J3(i-1,j-1) = P_calculated(i) - abs( V_abs(i)*V_abs(i)*abs(Y_matrix(i,i)) ) * cos( angle(Y_matrix(i,i)) );
            else
                J3(i-1,j-1) = -abs( V_abs(i)*V_abs(j)*abs(Y_matrix(i,j)) ) * cos( angle(Y_matrix(i,j)) + V_angle(j) - V_angle(i) );
            end
        end
    end

    [r1,c1]=size(J4);
    for i = 2 : r1+1
        for j = 2 : c1+1
            if i == j
                J4(i-1,j-1) = Q_calculated(i) - abs( V_abs(i)*V_abs(i)*abs(Y_matrix(i,i)) ) * sin( angle(Y_matrix(i,i)) );
            else
                J4(i-1,j-1) = -abs( V_abs(i)*V_abs(j)*abs(Y_matrix(i,j)) ) * sin( angle(Y_matrix(i,j)) + V_angle(j) - V_angle(i) );
            end
        end
    end
    
    J4_col=zeros(r1,1);
    J4_col(svc_bus-1)=-1;
    J4(:,(svc_bus-1))=J4_col;
    


    J=[J1 J2;J3 J4];

    % Correct Voltage Magnitudes and Angles
    Correction = inv(J) * Mismatch ;
    Correction=Correction';
    Correction(5:(4+svc_bus-2))= Correction(5:(4+svc_bus-2)).*V_abs(2:svc_bus-1);
    Correction((4+svc_bus):end)=Correction((4+svc_bus):end).*V_abs(svc_bus+1:end);

    V_abs(2:svc_bus-1) = V_abs(2:svc_bus-1) + Correction(5:(4+svc_bus-2));
    V_abs(svc_bus+1:end)=V_abs(svc_bus+1:end)+ Correction((4+svc_bus):end);

    V_angle(2:5) = V_angle(2:5) + Correction(1:4);
    Correction(4+svc_bus-1)=Correction(4+svc_bus-1);
    B_svc=B_svc+Correction(4+svc_bus-1);

  
    num_of_it=num_of_it+1;
    error=max(abs(Mismatch));
    % perror=abs(max(P_scheduled-P_calculated));
    % qerror=abs(max(Q_scheduled-Q_calculated));
    Correction
    B_svc
    Mismatch
    error
end
V_ANGLE=((V_angle)*180)/pi
V_abs'
disp(-V_abs(svc_bus)*V_abs(svc_bus)* B_svc*100*10e6)

Current_and_lineloss=zeros(7,4);
Current_and_lineloss(1,:)=current_and_lineloss(V_abs(1),V_ANGLE(1),V_abs(2),V_ANGLE(2),0.02,0.06);
Current_and_lineloss(2,:)=current_and_lineloss(V_abs(1),V_ANGLE(1),V_abs(3),V_ANGLE(3),0.08,0.24);
Current_and_lineloss(3,:)=current_and_lineloss(V_abs(2),V_ANGLE(2),V_abs(3),V_ANGLE(3),0.06,0.25);
Current_and_lineloss(4,:)=current_and_lineloss(V_abs(2),V_ANGLE(2),V_abs(4),V_ANGLE(4),0.06,0.18);
Current_and_lineloss(5,:)=current_and_lineloss(V_abs(2),V_ANGLE(2),V_abs(5),V_ANGLE(5),0.04,0.12);
Current_and_lineloss(6,:)=current_and_lineloss(V_abs(3),V_ANGLE(3),V_abs(4),V_ANGLE(4),0.01,0.03);
Current_and_lineloss(7,:)=current_and_lineloss(V_abs(4),V_ANGLE(4),V_abs(5),V_ANGLE(5),0.08,0.24);

fprintf('From To\t Current Mag.\tCurrent Angle(deg) Real Loss(MW) Reactive Loss(MVAr)\n');


% Display data
lines = {'1    2', '1    3', '2    3', '2    4', '2    5', '3    4', '4    5'};
for i = 1:size(Current_and_lineloss, 1)
    fprintf('%s\t\t%.4f\t\t%.4f\t\t%.6f\t\t%.6f\n', lines{i}, Current_and_lineloss(i, :));
end
