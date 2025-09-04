function S = build_scenarios(C, Nmc, varargin)
% Name-value args for convenience: build_scenarios(C,Nmc,'p_outage',0.05,'timeMode',"random")
p = inputParser;
addParameter(p,'p_outage',0.05);
addParameter(p,'timeMode',"random");
parse(p, varargin{:});
opts = p.Results;

% Load PV/Wind time series 
Tpv  = readtable(C.data.pv);    % expects Month, Hour, Irridiance/irradiance, Temp
Twd  = readtable(C.data.wind);  % expects Windspeed or Wind


pvIrrName = "Irridiance";
if ~ismember('Irridiance', Tpv.Properties.VariableNames)
    if ismember('Irradiance', Tpv.Properties.VariableNames)
        pvIrrName = "Irradiance";
    else
        error('PV table must have column Irridiance or Irradiance.');
    end
end

sigmaLoad = 0.05;  % normalized load noise

S = repmat(struct('month',[],'hour',[],'G',[],'TA',[],'Vw',[],'epsLoad',[],'outage',[]), Nmc,1);

for s=1:Nmc
    if opts.timeMode=="random"
        m = randi(12); h = randi([0 23]);
    else
        m = mod(s-1,12)+1; h = mod(s-1,24);
    end

    rows_pv = find(Tpv.Month==m & Tpv.Hour==h);
    if isempty(rows_pv), rows_pv = randi(height(Tpv)); end
    rpv = rows_pv(randi(numel(rows_pv)));
    G   = Tpv.(pvIrrName)(rpv);        % W/m^2
    TA  = Tpv.Temp(rpv);               % °C

    rwind = randi(height(Twd));
    if ismember('Windspeed',Twd.Properties.VariableNames)
        Vw = Twd.Windspeed(rwind);
    else
        Vw = Twd.Wind(rwind);
    end

    epsLoad = sigmaLoad*randn(1);

    if rand < opts.p_outage
        k   = randi(size(C.pairs,1));
        out = sprintf('%d-%d',C.pairs(k,1),C.pairs(k,2));
    else
        out = "none";
    end

    S(s) = struct('month',m,'hour',h,'G',G,'TA',TA,'Vw',Vw,'epsLoad',epsLoad,'outage',out);
end
end
