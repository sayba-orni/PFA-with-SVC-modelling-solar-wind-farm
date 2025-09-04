function Y5o = remove_line(Y5, C, a, b)
%REMOVE_LINE Remove the branch (a,b) from a 5-bus Ybus that was built
% with Y(a,b) = -y and per-end shunts jB at a and b.
% Inputs:
%   Y5  : 5x5 base Ybus (complex)
%   C   : struct with fields pairs, Zser, Bend
%   a,b : bus indices of the line to remove
% Output:
%   Y5o : modified Ybus with (a,b) branch removed

    % find the line (a,b) in the data
    idx = find( (C.pairs(:,1)==a & C.pairs(:,2)==b) | ...
                (C.pairs(:,1)==b & C.pairs(:,2)==a), 1);
    if isempty(idx)
        error('remove_line: line %d-%d not found in C.pairs.', a, b);
    end

    y = 1 / C.Zser(idx);   % series admittance
    B = C.Bend(idx);       % per-end shunt susceptance (jB) added at each end

    % start from the provided Y and subtract this branch's contributions
    Y5o = Y5;
    Y5o(a,a) = Y5o(a,a) - y - B;
    Y5o(b,b) = Y5o(b,b) - y - B;

    % off-diagonals were -y when the line was present → set to 0
    Y5o(a,b) = 0;
    Y5o(b,a) = 0;
end
