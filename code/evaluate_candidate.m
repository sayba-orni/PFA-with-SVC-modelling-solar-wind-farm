function R = evaluate_candidate(cand, C, S)
%EVALUATE_CANDIDATE  Run all scenarios for one SVC candidate (bus/mid-line).
% Produces per-scenario metrics: lossMW, voltage penalty, violation flag, etc.

    Y5_base = network5(C);
    Sbase   = C.Sbase;
    Vref    = C.Vref;
    Bmin    = C.Bmin;
    Bmax    = C.Bmax;

    Nsc     = numel(S);
    lossMW  = nan(Nsc,1);
    phiV    = nan(Nsc,1);
    vviol   = false(Nsc,1);
    qsvc    = nan(Nsc,1);
    conv    = false(Nsc,1);

    for s = 1:Nsc
        %  Outage-aware Ybus & edge list 
        if S(s).outage == "none"
            Y5 = Y5_base; out_ab = [];
            pairs_s = C.pairs; Zser_s = C.Zser; Bend_s = C.Bend;
        else
            ab      = sscanf(S(s).outage,'%d-%d'); 
            out_ab  = ab(:).';
            Y5      = remove_line(Y5_base, C, out_ab(1), out_ab(2));
            keep    = ~( (C.pairs(:,1)==out_ab(1) & C.pairs(:,2)==out_ab(2)) | ...
                         (C.pairs(:,1)==out_ab(2) & C.pairs(:,2)==out_ab(1)) );
            pairs_s = C.pairs(keep,:); 
            Zser_s  = C.Zser(keep); 
            Bend_s  = C.Bend(keep);
        end

        %  Islanding guard (base has 5 buses) 
        if is_islanded(pairs_s, 5)
            % Skip this scenario for all candidates (N/A, do not penalize)
            continue;
        end

        %  Scenario loads (pu) 
        Pd = (C.Pd_MW   .* (1 + S(s).epsLoad)) / Sbase;
        Qd = (C.Qd_MVAr .* (1 + S(s).epsLoad)) / Sbase;

        % Scheduled conventional gen (pu) 
        Pg = C.Pg_MW   / Sbase;
        Qg = C.Qg_MVAr / Sbase;

        % Add renewables: PV @ bus 5, Wind @ bus 4 (pu) 
        % (Assumes C.PV and C.WF are set in config.m)
        Ppv_MW   = pvPlantMW_datasheet(S(s).G,  S(s).TA,  C.PV);
        Pwind_MW = windPlantMW_formula(S(s).Vw, C.WF);
        Pg(5) = Pg(5) + Ppv_MW   / Sbase;   % bus 5 (PV)
        Pg(4) = Pg(4) + Pwind_MW / Sbase;   % bus 4 (Wind)

        %  Solve (bus or midline SVC) 
        if cand.type == "bus"
            res = solver_bus(Y5, Sbase, cand.k, Vref, C.SVC_Q_MAX, Pg, Qg, Pd, Qd);
        else
            % If the candidate line itself is outaged, skip (N/A)
            if ~isempty(out_ab) && isequal(sort(out_ab), sort([cand.a cand.b]))
                continue;
            end
            % Find that line in the current (post-outage) set
            eCase = find( (pairs_s(:,1)==cand.a & pairs_s(:,2)==cand.b) | ...
                          (pairs_s(:,1)==cand.b & pairs_s(:,2)==cand.a), 1);
            if isempty(eCase)
                continue;  % line not present after outage filter
            end
            res = solver_mid(pairs_s, Zser_s, Bend_s, Sbase, eCase, Vref, Bmin, Bmax, Pg, Qg, Pd, Qd);
        end

        % Convergence flag 
        conv(s) = isfield(res,'converged') && res.converged;

        %  reject numerically-wild solutions (treat like infeasible) -
        Vabs = res.V_abs(:);
        if ~conv(s) || any(~isfinite(Vabs)) || any(Vabs < 0.8 | Vabs > 1.2)
            lossMW(s) = inf;
            phiV(s)   = inf;
            vviol(s)  = true;
            qsvc(s)   = NaN;
            conv(s)   = false;
            continue;
        end

        % If we're here, it's sane; record metrics
        lossMW(s) = res.totPloss;

        % voltage penalty on buses 2..end (exclude slack only)
        if numel(Vabs) >= 2
            Vcheck = Vabs(2:end);
        else
            Vcheck = Vabs;  % defensive, shouldn't happen
        end
        [phiV(s), vviol(s)] = volt_penalty(Vcheck, C.Vband);

        qsvc(s) = res.Qsvc_MVAr;
    end

    R = struct('lossMW',lossMW, 'phiV',phiV, 'vviol',vviol, ...
               'qsvc',qsvc, 'converged',conv);
end
