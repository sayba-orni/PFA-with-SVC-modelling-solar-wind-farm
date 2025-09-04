function Y5 = network5(C)
Y5 = zeros(5,5);
for e=1:size(C.pairs,1)
  a=C.pairs(e,1); b=C.pairs(e,2); y=1/C.Zser(e); B=C.Bend(e);
  Y5(a,a)=Y5(a,a)+y+B; Y5(b,b)=Y5(b,b)+y+B; Y5(a,b)=Y5(a,b)-y; Y5(b,a)=Y5(b,a)-y;
end
end

function Y5o = remove_line(Y5, C, a, b)
% remove series y and end shunts for line (a,b)
idx = find( (C.pairs(:,1)==a & C.pairs(:,2)==b) | (C.pairs(:,1)==b & C.pairs(:,2)==a), 1);
y = 1/C.Zser(idx); B = C.Bend(idx);
Y5o = Y5;
Y5o(a,a)=Y5o(a,a)-y-B; Y5o(b,b)=Y5o(b,b)-y-B; Y5o(a,b)=0; Y5o(b,a)=0;
end
