function Zpv = GenZpv(pv,nb,nl,f,Zbranch)
%   算PV的回路阻抗矩阵
npv = length(pv); % PV节点数
Zpv = zeros(npv,1); % PV阻抗矩阵 npv*1大小
for ipv = 1:npv
    V = zeros(nb,1);
    I = zeros(nb,1);
    I(pv(ipv)) = -1; % 在其中一条PV节点上注入-1A的电流时，阻抗等于其总线电压。
    for k = nl:-1:2 % 回代过程
        i = f(k);
        I(i,1) = I(i,1) + I(k,1);
    end
    for k = 2:1:nl    % 前推过程
        i = f(k);
        V(k,1) = V(i,1) - Zbranch(k,1) * I(k,1);
    end 
    Zpv(ipv,1) = V(pv(ipv));
end

