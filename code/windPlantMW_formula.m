function Pmw = windPlantMW_formula(v, W)
% Piecewise curve with a*v^3 + b*Pr between vin and vr, Pr at vr..vout
a = W.PrMW / (W.v_r^3 - W.v_in^3);
b = (W.v_in^3) / (W.v_r^3 - W.v_in^3);

Pmw = zeros(size(v));
Pmw(v < W.v_in | v > W.v_out) = 0;
mask = (v >= W.v_in) & (v <= W.v_r);
Pmw(mask) = a * v(mask).^3 + b * W.PrMW;
Pmw(v > W.v_r & v <= W.v_out) = W.PrMW;
end
