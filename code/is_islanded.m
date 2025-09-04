function tf = is_islanded(pairs, Nnodes)
%IS_ISLANDED  True if the network defined by "pairs" is disconnected.
% pairs  : Ex2 edge list over integer node labels 1..Nnodes
% Nnodes : total number of buses in the base system (here, 5)

    if nargin < 2 || isempty(Nnodes)
        % Best-effort guess (not recommended for our case)
        Nnodes = max(pairs(:));
    end

    if isempty(pairs)
        tf = true;  % no edges at all ⇒ disconnected
        return;
    end

    % Important: tell graph the total node count so isolated nodes are included.
    G = graph(pairs(:,1), pairs(:,2), [], Nnodes);
    comps = conncomp(G);      % component id for each node
    tf = (max(comps) > 1);    % more than one component ⇒ islanded
end
