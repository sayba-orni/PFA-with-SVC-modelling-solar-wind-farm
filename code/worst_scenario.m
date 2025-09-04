%% ------------------------------------------------------------
% Compare voltage profiles: MID-4-5 vs BUS-5 in worst scenarios (robust)
% ------------------------------------------------------------
candA = struct('type',"mid",'k',[],'a',4,'b',5);   % MID-4-5
candB = struct('type',"bus",'k',5,'a',[],'b',[]);  % BUS-5

% Evaluate all scenarios once
RA = evaluate_candidate(candA, C, S);
RB = evaluate_candidate(candB, C, S);

% Helper: index of worst VALID scenario (converged & finite)
worst_valid_idx = @(R,lambda) ...
    (function ()
        valid = R.converged & isfinite(R.lossMW) & isfinite(R.phiV);
        if any(valid)
            Jv = R.lossMW(valid) + lambda*R.phiV(valid);
            [~,k] = max(Jv);
            idxs = find(valid);
            idx = idxs(k);
        else
            idx = NaN;
        end
    end)();

idxA_worst = worst_valid_idx(RA, lambda);
idxB_worst = worst_valid_idx(RB, lambda);

% Log what we picked
fprintf('\n-- Worst VALID scenario indices --\n');
if ~isnan(idxA_worst)
    JA = RA.lossMW + lambda*RA.phiV;
    fprintf('MID-4-5 worst s = %d, score %.4f (loss %.4f, phi %.6f)\n', ...
        idxA_worst, JA(idxA_worst), RA.lossMW(idxA_worst), RA.phiV(idxA_worst));
else
    fprintf('MID-4-5 has no valid scenarios (all skipped / non-converged).\n');
end
if ~isnan(idxB_worst)
    JB = RB.lossMW + lambda*RB.phiV;
    fprintf('BUS-5   worst s = %d, score %.4f (loss %.4f, phi %.6f)\n', ...
        idxB_worst, JB(idxB_worst), RB.lossMW(idxB_worst), RB.phiV(idxB_worst));
else
    fprintf('BUS-5 has no valid scenarios (very unlikely here).\n');
end

run_one = @(cand, sidx) run_one_scenario(cand, C, S(sidx));

% Own-worst runs (only if index is valid)
resA_worstA = struct('converged',false);
resB_worstB = struct('converged',false);
if ~isnan(idxA_worst), resA_worstA = run_one(candA, idxA_worst); end
if ~isnan(idxB_worst), resB_worstB = run_one(candB, idxB_worst); end

% Apples-to-apples: use BUS-5's worst valid scenario if present,
% otherwise use MID-4-5's worst valid scenario
common_idx = idxB_worst;
if isnan(common_idx), common_idx = idxA_worst; end
if isnan(common_idx)
    error('No common valid scenario to compare.');
end

resA_onCommon = run_one(candA, common_idx);
resB_onCommon = run_one(candB, common_idx);

% Safe extractor for |V| at buses 2..5
takeV2to5 = @(res) ( ...
    (isstruct(res) && isfield(res,'V_abs') && ~isempty(res.V_abs)) ...
    * ones(1,0) );  % dummy to allow function handle

function V = takeV2to5_impl(res)
    if ~isstruct(res) || ~isfield(res,'V_abs') || isempty(res.V_abs)
        V = nan(1,4);
        return
    end
    V = res.V_abs(2:min(5,numel(res.V_abs)));
    if numel(V) < 4, V(end+1:4) = NaN; end
end
% wrap to call the local subfunction above
takeV2to5 = @(res) takeV2to5_impl(res);

V_A_own   = takeV2to5(resA_worstA);
V_B_own   = takeV2to5(resB_worstB);
V_A_same  = takeV2to5(resA_onCommon);
V_B_same  = takeV2to5(resB_onCommon);
busTicks  = 2:5; busLabels = compose('Bus %d', busTicks);

if ~exist('figures','dir'), mkdir figures; end

% Figure 1: each on its own worst (plot only if at least one valid)
if any(isfinite(V_A_own)) || any(isfinite(V_B_own))
    f1 = figure('Name','Voltage profile: each on its own worst scenario', ...
                'Color','w','NumberTitle','off','Visible','on');
    hold on
    if any(isfinite(V_A_own)), plot(busTicks, V_A_own, '-o','LineWidth',1.5); end
    if any(isfinite(V_B_own)), plot(busTicks, V_B_own, '-s','LineWidth',1.5); end
    yline(C.Vband(1),'--'); yline(C.Vband(2),'--');
    grid on; xlim([2 5]); ylim([0.94 1.06]);
    xticks(busTicks); xticklabels(busLabels);
    ylabel('|V| (pu)'); title('Own-worst scenarios');
    legendStr = {'MID-4-5 (own worst)','BUS-5 (own worst)','Vmin','Vmax'};
    legend(legendStr{1:2+(~isempty(C.Vband))*2},'Location','best'); %#ok<NASGU>
    exportgraphics(f1, fullfile('figures','volts_own_worst.png'), 'Resolution',200);
end

% Figure 2: both on the same (common) worst scenario
f2 = figure('Name','Voltage profile: both on common worst scenario', ...
            'Color','w','NumberTitle','off','Visible','on');
hold on
plot(busTicks, V_A_same, '-o','LineWidth',1.5);
plot(busTicks, V_B_same, '-s','LineWidth',1.5);
yline(C.Vband(1),'--'); yline(C.Vband(2),'--');
grid on; xlim([2 5]); ylim([0.94 1.06]);
xticks(busTicks); xticklabels(busLabels);
ylabel('|V| (pu)'); title(sprintf('Same scenario (s=%d)', common_idx));
legend('MID-4-5','BUS-5','Vmin','Vmax','Location','best');
exportgraphics(f2, fullfile('figures','volts_same_scenario.png'), 'Resolution',200);

% Print summary for the common scenario
fprintf('\n-- Detailed metrics on common scenario (s=%d) --\n', common_idx);
fprintf('MID-4-5: converged=%d, loss=%.4f MW, Qsvc=%.2f MVAr, mode=%s\n', ...
    isfield(resA_onCommon,'converged') && resA_onCommon.converged, ...
    getfield(resA_onCommon,'totPloss',nan), getfield(resA_onCommon,'Qsvc_MVAr',nan), ...
    getfield(resA_onCommon,'mode','N/A'));
fprintf('BUS-5  : converged=%d, loss=%.4f MW, Qsvc=%.2f MVAr, mode=%s\n', ...
    isfield(resB_onCommon,'converged') && resB_onCommon.converged, ...
    getfield(resB_onCommon,'totPloss',nan), getfield(resB_onCommon,'Qsvc_MVAr',nan), ...
    getfield(resB_onCommon,'mode','N/A'));
