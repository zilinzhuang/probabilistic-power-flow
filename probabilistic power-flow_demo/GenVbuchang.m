function V_compen = GenVbuchang(I_compen,pv,nb,nl,f,Zb)
%  含PV节点时的电压矫正
V_compen = zeros(nb,1); % 
I = zeros(nb,1);
I(pv) = I_compen;
% 回代过程
for k = nl:-1:2
    i = f(k);
    I(i) = I(i) + I(k);
end
% 前推过程
for k = 2:nl
    i = f(k);
    V_compen(k) = V_compen(i) - Zb(k) * I(k);
end
