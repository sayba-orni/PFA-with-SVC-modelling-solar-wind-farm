function Pmw = pvPlantMW_datasheet(G_Wm2, TA_C, P)
% P – struct from C.PV
% Implements:
%  FF = (VMPP*IMPP)/(VOC*ISC)
%  Tcy = TA + s*((NOT-20)/0.8)
%  Vy  = VOC - Kv*Tcy
%  Iy  = s * ( ISC + Ki*(Tcy-25) )
%  Po  = N * FF * Vy * Iy
s = max(0, G_Wm2/1000);                 % kW/m^2 (0..~1), night => 0
FF = (P.VMPP*P.IMPP)/(P.VOC*P.ISC);
Tcy = TA_C + s*((P.NOT-20)/0.8);
Vy  = P.VOC - P.Kv*Tcy;
Iy  = s * ( P.ISC + P.Ki*(Tcy-25) );
Po_W = P.N_modules * FF * Vy .* Iy;     % Watts
Pmw  = max(0, Po_W) / 1e6;              % MW, clamp negative to 0
end
