clc; clear; close all;
data33; % 负荷、支路阻抗数据均为标幺值
SB = 10; % MVA 基准功率
SB = SB*1000; % 基准功率转换为kVA
UB = 12.66;   % 基准电压 kV
IB = SB/(sqrt(3)*UB); % 基准电流 A

U_ref = 1.0; % 首节点/平衡节点
NB = 33; NL = 32; % 总节点数, 总线路数
fr = mpc.branch(:,1); % 首节点行向量 from
Zbranch = mpc.branch(:,3)+1i*mpc.branch(:,4); % p.u.
Ybranch = 1./Zbranch; % p.u.
S_d = mpc.bus(:,3)+1i*mpc.bus(:,4); % 负荷 p.u.
Sd_orginal = S_d; % 原负荷数据
ep = 1e-2; iter_max = 100; % 收敛精度,最大迭代次数
U = U_ref * ones(NB,1);
Uold = U;
iter_num = 0;

fr = [0; fr]; % 添加一条虚拟支路至所有支路的顶层, 这样支路数 就等于 收节点数，以便用于后面迭代计算
Zbranch = [0; Zbranch];
NL = NL + 1;
% PQ型 DG
Ind_pq = 9; % PQ型DG接入位置
Ppq_G = 300/SB; % kW -> p.u.
Qpq_G = 100/SB; % kVar -> p.u
% PI型 DG
Ind_pi = 20; % PI型DG接入位置
Ppi_G = 100/SB; % kW -> p.u.
Ipi_G = 50/IB; % A -> p.u
% PV型 DG
Ind_pv = 32; % PV型DG 接入位置
Ppv_G = 600/SB; % kW -> p.u
Qpv_G = 600/SB; % kVar PV型DG的无功出力调节范围 -> p.u
Upv = [1.0]; % PV型DG 定电压
Ug = Upv(length(Ind_pv));

S_d(Ind_pq) = S_d(Ind_pq) - complex(Ppq_G,Qpq_G);% 对 PQ型DG 将其作为负的负荷添加到节点负荷处
S_d(Ind_pi) = S_d(Ind_pi) - Ppi_G;               % 对 PI型DG 先将其恒定的有功输出作为负的负荷添加到节点负荷处
S_d(Ind_pv) = S_d(Ind_pv) - Ppv_G;               % 对 PV型DG 先将其恒定的有功输出作为负的负荷添加到节点负荷处
% 生成PV节点的阻抗矩阵
Zpv = GenZpv(Ind_pv,NB,NL,fr,Zbranch); % 子函数
if size(Zpv,1) > 0
	Bpv = (imag(Zpv))^-1;
end
% 电压校正潮流
npv = length(Ind_pv);
Qpv = zeros(npv,1);
while iter_num < iter_max
	iter_num = iter_num + 1;
	S = S_d;
	S_to = S;
	S_fr = S_to;
	% 回代过程
	for k = NL:-1:2 % 从最后一条支路,依次-1, 至第一条支路
		i = fr(k);
		S_fr(k) = S_to(k) + Zbranch(k) * abs(S_to(k)/U(k))^2;
		S_to(i) = S_to(i) + S_fr(k);
	end
	% 前推
	for k = 2:NL % 从第一条支路,依次+1, 至最后一条支路
		i = fr(k);
		U(k) = U(i) - Zbranch(k) * conj(S_fr(k)/U(i));
	end
	% 检查收敛
	dU = abs(U - Uold); % 偏差量
	
	Uold = U;
	if ~isempty(Ind_pv) % 存在PV节点
		% 计算PV节点的无功功率校正
		dE = ( Ug./abs(U(Ind_pv))-1).*real(U(Ind_pv) ); %
		dD = Bpv * dE;
		dC = dD .* imag(U(Ind_pv))./real(U(Ind_pv));
		V_corr = GenVcorr(dC+1j*dD, Ind_pv, NB, NL, fr, Zbranch); % 生成节点电压矫正量 子函数
		U = U + V_corr;                          % 代入节点电压迭代中
		DQ = dD .* abs(U(Ind_pv)).^2 ./ real(U(Ind_pv)); % PV节点 无功输出量
		Qpv = Qpv + DQ;                          % 更新 PV节点 无功功率输出
		S_d(Ind_pv) = S_d(Ind_pv) - 1j*DQ;               % 将PV节点输出的无功 代入节点功率中用于计算
	end
	
	if ~isempty(Ind_pi) % 存在PI节点
		Q_PI = sqrt( Ipi_G^2 * abs(real(U(Ind_pi))^2 + imag(U(Ind_pi))^2) - Ppi_G^2);
		if abs(Q_PI) > Ppi_G % 判断无功输出是否越界，一般无功输出不能超过有功输出范围，可视为有功输出最大值的正负范围内，即 QG = [-PG,+PG]
			Q_PI = Ppi_G;
		end
		S_d(Ind_pi) = S_d(Ind_pi) - 1j*Q_PI;  % 将PI节点计算得出的无功 代入节点功率中用于计算
		
	end
	if max(dU) <= ep % 未达到收敛精度
		break
	end
end
% 保存发电机无功输出 和 支路功率流
if ~isempty(Ind_pv)
	Qpv_G = Qpv; % p.u.
end
Sslack = S_to(1);
S_fr = S_fr(2:end); % i->j功率
S_to = S_to(2:end); % j->i功率
Sloss = S_fr - S_to;
Ploss = SB*real(Sloss); % kW 有功网损
Qloss = SB*imag(Sloss); % kVar 无功网损


