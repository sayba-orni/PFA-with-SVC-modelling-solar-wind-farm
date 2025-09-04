function T = rank_candidates(C, scenarios, candidates, lambda, alpha)
% lambda: weight on voltage penalty (MW per pu-sum), alpha: CVaR tail
M = numel(candidates);

meanLoss  = nan(M,1);
meanPhi   = nan(M,1);
violProb  = nan(M,1);
cvarScore = nan(M,1);
meanScore = nan(M,1);

bigM = 5e3;      % fallback score if a candidate has 0 valid scenarios
minWarnFrac = 0.50;  % warn if fewer than 50% scenarios are valid

for m = 1:M
    % name (for logs and table)
    if candidates(m).type == "bus"
        nm = sprintf('BUS-%d', candidates(m).k);
    else
        nm = sprintf('MID-%d-%d', candidates(m).a, candidates(m).b);
    end

    % evaluate all scenarios for this candidate 
    R = evaluate_candidate(candidates(m), C, scenarios);

    % valid scenarios = converged & finite metrics
    valid = R.converged & isfinite(R.lossMW) & isfinite(R.phiV);

    % log feasibility
    nSc   = numel(valid);
    nBad  = nSc - nnz(valid);
    fprintf('%-10s  infeasible/nonconv: %3d / %3d (%.1f%%)\n', ...
            nm, nBad, nSc, 100*nBad/max(1,nSc));

    if ~any(valid)
        % no usable scenarios → assign fallback scores
        meanLoss(m)  = bigM;
        meanPhi(m)   = bigM;
        violProb(m)  = 1.0;
        cvarScore(m) = bigM;
        meanScore(m) = bigM;
        continue;
    end

    % Optional: warning if too few valid scenarios
    if (nnz(valid) < minWarnFrac*nSc)
        warning('%s: only %d/%d (%.0f%%) scenarios valid; rankings may be noisy.', ...
                 nm, nnz(valid), nSc, 100*nnz(valid)/nSc);
    end

    %  metrics on valid scenarios only 
    loss  = R.lossMW(valid);
    phi   = R.phiV(valid);
    viol  = R.vviol(valid);

    % combined objective per scenario
    J = loss + lambda*phi;

    % defensive (should be finite already)
    loss(~isfinite(loss)) = [];
    phi(~isfinite(phi))   = [];
    J(~isfinite(J))       = [];

    % aggregate stats
    meanLoss(m)  = mean(loss);
    meanPhi(m)   = mean(phi);
    violProb(m)  = mean(double(viol));
    meanScore(m) = mean(J);

    % CVaR of the valid objective values
    cvarScore(m) = cvar(J, alpha);
end

% build output table 
name = strings(M,1);
for m = 1:M
    if candidates(m).type == "bus"
        name(m) = sprintf('BUS-%d', candidates(m).k);
    else
        name(m) = sprintf('MID-%d-%d', candidates(m).a, candidates(m).b);
    end
end

T = table(name, meanLoss, meanPhi, violProb, meanScore, cvarScore);
T = sortrows(T, 'cvarScore');
end
