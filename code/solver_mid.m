function res = solver_mid(pairs, Zser, Bend, Sbase, eCase, Vref, Bmin, Bmax, Pg, Qg, Pd, Qd)
% SOLVER_MID  Robust NR power flow with an SVC at the MIDPOINT of line eCase.
% Splits the selected line into two half-π sections with new bus 6 (|V6|=Vref).
%
% Inputs
%   pairs, Zser, Bend : line endpoints, series impedances, and per-end shunts
%   Sbase             : base MVA
%   eCase             : index of the line to split (row in 'pairs')
%   Vref              : |V6| setpoint
%   Bmin,Bmax         : SVC susceptance bounds (pu)
%   Pg,Qg,Pd,Qd       : 1x5 vectors (pu) for the original buses (1..5)
%
% Output (same fields as solver_bus, but with N=6 and bus 6 as SVC)

    % build base 5-bus Ybus 
    Y5 = zeros(5,5);
    E  = size(pairs,1);
    for e = 1:E
        a = pairs(e,1); b = pairs(e,2);
        y = 1/Zser(e); B = Bend(e);
        Y5(a,a)=Y5(a,a)+y+B; Y5(b,b)=Y5(b,b)+y+B;
        Y5(a,b)=Y5(a,b)-y;   Y5(b,a)=Y5(b,a)-y;
    end

    % select line to split 
    a = pairs(eCase,1); b = pairs(eCase,2);
    Zab = Zser(eCase);  Bab = Bend(eCase);
    y_ab = 1/Zab;

    % create 6-bus Y by splitting line (a,b)
    Y = zeros(6,6); Y(1:5,1:5) = Y5;

    % remove original (a,b) branch (series and its end shunts)
    Y(a,a) = Y(a,a) - y_ab - Bab;
    Y(b,b) = Y(b,b) - y_ab - Bab;
    Y(a,b) = 0;  Y(b,a) = 0;

    % add two half-π sections: a–6 and 6–b
    z_half = Zab/2;   y_half = 1/z_half;   B_half = Bab/2;
    Y(a,a) = Y(a,a) + y_half + B_half;
    Y(6,6) = Y(6,6) + y_half + B_half;
    Y(a,6) = Y(a,6) - y_half;  Y(6,a) = Y(6,a) - y_half;

    Y(b,b) = Y(b,b) + y_half + B_half;
    Y(6,6) = Y(6,6) + y_half + B_half;
    Y(b,6) = Y(b,6) - y_half;  Y(6,b) = Y(6,b) - y_half;

    % NR setup
    N=6; m=6;
    Vmag = ones(1,N); Vmag(1)=1.06; Vmag(2:5)=1.00; Vmag(m)=Vref;
    Vang = zeros(1,N);
    Bsvc = 0.0; mode="REG";

    % schedules (append 0 at bus 6)
    Psch = [Pg 0] - [Pd 0];
    Qsch = [Qg 0] - [Qd 0];

    tol=1e-6; max_iter=50; iter=0; err=1;

    while err>tol && iter<max_iter
        iter=iter+1;

        Yit = Y; Yit(m,m) = Yit(m,m) + 1i*Bsvc;

        [P,Q] = pq_injections(Vmag,Vang,Yit);
        dP = Psch - P; dQ = Qsch - Q;
        M  = [dP(2:N) dQ(2:N)]';
        err = max(abs(M)); if err<=tol, break; end

        [J1,J2,J3,J4] = jac_blocks(Vmag,Vang,Yit);
        if mode=="REG"
            idx=m-1;
            J2(:,idx)=0; J4(:,idx)=0; J4(idx,idx) = -Vmag(m)^2;
        end
        J=[J1 J2; J3 J4];

        dx = (J + 1e-12*eye(size(J))) \ M; dx=dx.';
        PQ=N-1; dth=dx(1:PQ); dVm=dx(PQ+1:end);
        dB = (mode=="REG") * dVm(m-1);

        % backtracking + bounds
        alpha=1.0; accepted=false; M0=0.5*(M.'*M);
        while alpha>1/1024
            Vang_try=Vang;  Vang_try(2:N)=Vang(2:N)+alpha*dth;
            dV_full=zeros(1,N); dV_full(2:N)=alpha*dVm;
            Vmag_try=Vmag + dV_full; if mode=="REG", Vmag_try(m)=Vref; end

            B_try=Bsvc + alpha*dB; mode_try=mode;
            if mode=="REG" && (B_try<Bmin || B_try>Bmax)
                B_try=min(max(B_try,Bmin),Bmax); mode_try="SAT";
            end

            Ytr=Y; Ytr(m,m)=Ytr(m,m)+1i*B_try;
            [Pt,Qt]=pq_injections(Vmag_try,Vang_try,Ytr);
            Mt=[ (Psch(2:N)-Pt(2:N))' ; (Qsch(2:N)-Qt(2:N))' ];
            if 0.5*(Mt.'*Mt) < M0
                Vang=Vang_try; Vmag=Vmag_try; Bsvc=B_try; mode=mode_try; accepted=true; break;
            else
                alpha=alpha/2;
            end
        end

        if ~accepted
            Vang(2:N)=Vang(2:N)+1e-3*dth;
            if mode=="REG"
                Vmag(2:N)=Vmag(2:N)+1e-3*dVm; Vmag(m)=Vref;
                Bsvc=min(max(Bsvc + 1e-3*dB, Bmin), Bmax);
            else
                Vmag(2:N)=Vmag(2:N)+1e-3*dVm;
            end
        end

        if mode=="SAT"
            if Bsvc>Bmin+1e-6 && Bsvc<Bmax-1e-6 && abs(Vmag(m)-Vref)<5e-4
                mode="REG"; Vmag(m)=Vref;
            end
        end
    end

    % post-processing 
    Yfinal=Y; Yfinal(m,m)=Yfinal(m,m)+1i*Bsvc;
    Vdeg=Vang*180/pi;
    Qsvc_MVAr = -(Vmag(m)^2)*Bsvc*Sbase; svcType='Inductive'; if Qsvc_MVAr>=0, svcType='Capacitive'; end

    % branch list for currents/losses (replace a-b by a-6 and 6-b)
    basePairs=pairs; baseZ=Zser;
    keep = ~( (basePairs(:,1)==a & basePairs(:,2)==b) | (basePairs(:,1)==b & basePairs(:,2)==a) );
    basePairs = basePairs(keep,:); baseZ = baseZ(keep);
    lines = [basePairs; a 6; 6 b];
    Zs    = [baseZ;     z_half; z_half];

    L=size(lines,1);
    I_mag=zeros(L,1); I_ang=zeros(L,1); Pl=zeros(L,1); Ql=zeros(L,1);
    for ee=1:L
        ii=lines(ee,1); kk=lines(ee,2);
        [I_mag(ee),I_ang(ee),Pl(ee),Ql(ee)] = current_and_lineloss(Vmag(ii),Vdeg(ii),Vmag(kk),Vdeg(kk),Zs(ee),Sbase);
    end

    % slack injections
    Pg1_MW = Sbase*calc_injection_P(1,Vmag,Vang,Yfinal);
    Qg1_MVAr = Sbase*calc_injection_Q(1,Vmag,Vang,Yfinal);

    % outputs
    res.iter=iter; res.err=err; res.converged=(err<=tol); res.mode=char(mode);
    res.svc_bus=m; res.B_svc=Bsvc; res.Qsvc_MVAr=Qsvc_MVAr; res.svcType=svcType;
    res.V_abs=Vmag; res.V_deg=Vdeg;

    res.genMW    = [Pg1_MW,   Pg(2:end)*Sbase, 0];
    res.genMVAr  = [Qg1_MVAr, Qg(2:end)*Sbase, 0];
    res.loadMW   = [Pd 0] * Sbase;
    res.loadMVAr = [Qd 0] * Sbase;

    res.lines=lines; res.I_mag=I_mag; res.I_ang=I_ang;
    res.Lmw=Pl; res.Lmvar=Ql; res.totPloss=sum(Pl); res.totQloss=sum(Ql);
end

% helpers
function [P,Q] = pq_injections(Vmag, Vang, Y)
    N=numel(Vmag); P=zeros(1,N); Q=zeros(1,N);
    for i=1:N
        for k=1:N
            th = angle(Y(i,k)) + Vang(k) - Vang(i);
            g  = abs(Y(i,k));
            P(i) = P(i) + Vmag(i)*Vmag(k)*g*cos(th);
            Q(i) = Q(i) - Vmag(i)*Vmag(k)*g*sin(th);
        end
    end
end

function [J1,J2,J3,J4] = jac_blocks(Vmag,Vang,Y)
    N=numel(Vmag); n=N-1;
    J1=zeros(n,n); J2=zeros(n,n); J3=zeros(n,n); J4=zeros(n,n);
    [P,Q] = pq_injections(Vmag,Vang,Y);
    for i=2:N
        for k=2:N
            if i==k
                Yii = Y(i,i); gii=abs(Yii); thi=angle(Yii);
                J1(i-1,k-1) = -Q(i) - Vmag(i)^2*gii*sin(thi);
                J2(i-1,k-1) =  P(i) + Vmag(i)^2*gii*cos(thi);
                J3(i-1,k-1) =  P(i) - Vmag(i)^2*gii*cos(thi);
                J4(i-1,k-1) =  Q(i) - Vmag(i)^2*gii*sin(thi);
            else
                th = angle(Y(i,k)) + Vang(k) - Vang(i);
                g  = abs(Y(i,k));
                J1(i-1,k-1) = -Vmag(i)*Vmag(k)*g*sin(th);
                J2(i-1,k-1) =  Vmag(i)*Vmag(k)*g*cos(th);
                J3(i-1,k-1) = -Vmag(i)*Vmag(k)*g*cos(th);
                J4(i-1,k-1) = -Vmag(i)*Vmag(k)*g*sin(th);
            end
        end
    end
end

function [Iabs,Iang,Ploss,Qloss] = current_and_lineloss(Vm_i,ang_i_deg,Vm_j,ang_j_deg,Z,Sbase)
    Vi = Vm_i * exp(1i*deg2rad(ang_i_deg));
    Vj = Vm_j * exp(1i*deg2rad(ang_j_deg));
    I  = (Vi - Vj) / Z;
    Iabs = abs(I); Iang = rad2deg(angle(I));
    S_loss_pu = (Iabs^2) * Z;      % I^2 * Z
    Ploss = Sbase * real(S_loss_pu);
    Qloss = Sbase * imag(S_loss_pu);
end

function P = calc_injection_P(i,Vmag,Vang,Y)
    N=numel(Vmag); P=0;
    for k=1:N
        th = angle(Y(i,k)) + Vang(k) - Vang(i);
        P = P + Vmag(i)*Vmag(k)*abs(Y(i,k))*cos(th);
    end
end

function Q = calc_injection_Q(i,Vmag,Vang,Y)
    N=numel(Vmag); Q=0;
    for k=1:N
        th = angle(Y(i,k)) + Vang(k) - Vang(i);
        Q  = Q - Vmag(i)*Vmag(k)*abs(Y(i,k))*sin(th);
    end
end
