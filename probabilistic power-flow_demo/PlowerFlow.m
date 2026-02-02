clc
clear
close all

%% 假设条件

MODEL = data33; % 负荷、支路阻抗数据均为 标幺值

SB = 10;        % MVA 基准功率
SB = SB*1000;   % 基准功率转换为kVA
UB = 12.66;     % 基准电压 kV
IB = SB/(sqrt(3)*UB); % 基准电流 A


U_ref = 1.0; % 首节点/平衡节点
NUMBUS = 33; NUMLINE = 32; % 总节点数, 总线路数
firstbranchbus = MODEL.branch(:,1); % 首节点列向量
Zbranch = MODEL.branch(:,3)+1i*MODEL.branch(:,4); % 支路阻抗矩阵p.u.
Ybranch = 1./Zbranch;                         % 支路的电纳矩阵p.u.
loadbus = MODEL.bus(:,3)+1i*MODEL.bus(:,4);           % 负荷阻抗矩阵 p.u.
Sd_orginal = loadbus;                             % 原负荷数据
ep = 1e-5; 
iter_max = 100;                 % 收敛精度及最大迭代次数？
U = U_ref * ones(NUMBUS,1);
Uold = U;
iter_num = 0;



firstbranchbus = [0; firstbranchbus]; % 添加一条虚拟支路至所有支路的顶层, 这样支路数就等于节点数，便用于后面迭代计算
Zbranch = [0; Zbranch];    %假设与首节点连接的支路阻抗为0，并把支路数加1.
NUMLINE = NUMLINE + 1;             
Uppf=zeros(NUMBUS,1);                     %概率潮流电压
Ucrossline=zeros(NUMBUS,1);                % 电压越限数

%% Monte-Carlo概率模型建立

NMC=2000;                          %  monte-carlo模拟数

%% LHS
windLHS=wblinv([0 : 1/NMC : 1],7,2);      %LHS抽样开始！
sunLHS=betainv([0 : 1/NMC : 1],4,0.1);

Index1 = randperm(NMC);   % 生成随机排列数1
Index2 = randperm(NMC);   % 生成随机排列数2
windLHS = windLHS(Index1);      % 随机排列后的矩阵
sunLHS = sunLHS(Index2);      % 随机排列后的矩阵
pw1=wind2power(windLHS,8,1)*1000;       %  风能有功功率
ps1=solar2power(sunLHS,1000)*1000;       %  太阳能功率

%% 随机抽样MC
%%wind1=wblrnd(7,2,1,NMC);                %  设定WEIBULL参数！scale parameter A and shape parameter B.后两项为矩阵形式
%%sun1=betarnd(4,0.1,1,NMC);             %  设定BETA参数！
%%pw2=wind2power(wind1,5,1)*1000;       %  风能有功功率
%%ps2=solar2power(sun1,1000)*1000;       %  太阳能功率




fprintf('请输入需要哪个节点的电压概率图像？(范围1-33)\n');
printnumbus=input(['节点号=']);
Uprintbus=zeros(NMC,NUMBUS);          %建立所有模拟的数据库矩阵

plossppf=zeros(1,NMC);                %有功网损
RANDOMMC = SUIJIMCF(printnumbus,NMC); %随机抽样MC
LHSMCMC =  LHSMC(printnumbus,NMC);    %LHSMC
time0 = clock;  % 计时开始

for mci=1:NMC 
%初始化
loadbus = Sd_orginal;    % 负荷阻抗矩阵 p.u.
U = U_ref * ones(NUMBUS,1);
Uold = U;
iter_num = 0; 

%% 对DG分类、定义

% PQ型 DG        风机
PQNUM = [33,15]; % PQ型DG位置
Ppq_G = pw1(1,mci)/SB/2; % kW >>> p.u.
Qpq_G = Ppq_G.*tan(acos(0.98)); % kVar >>> p.u  如何更精确的计算呢 and 风机分类
% PI型 DG
PINUM = [];    % PI型DG位置
Ppi_G = 20/SB; % kW >>> p.u.
Ipi_G = -10/IB; % A >>> p.u
% PV型 DG         光伏、燃料电池
PVNUM = [18,8]; % PV型DG位置
Ppv_G = ps1(1,mci)/SB/2; % kW >>> p.u
Qpv_G = 0/SB; % PV型DG的无功出力调节范围 kVar>>> p.u
Upv = [0.98,0.98]; % PV型DG 定电压


loadbus(PQNUM) = loadbus(PQNUM) - complex(Ppq_G,Qpq_G);% 对PQ DG 把其作为负的负荷加到节点负荷处
loadbus(PINUM) = loadbus(PINUM) - Ppi_G;               % 对PI DG 先将其恒定的有功输出作为负的负荷添加到节点负荷处
loadbus(PVNUM) = loadbus(PVNUM) - Ppv_G;               % 对PV DG 先将其恒定的有功输出作为负的负荷添加到节点负荷处


% 生成PV节点的阻抗矩阵 
Zpv = GenZpv(PVNUM,NUMBUS,NUMLINE,firstbranchbus,Zbranch); % 子函数
if (size(Zpv,1) > 0)
	Bpv = (imag(Zpv))^-1;              % 虚部电纳
end




%% 前推回代潮流
npv = length(PVNUM);
Qpv = zeros(npv,1);                   %PV节点无功

while iter_num < iter_max             %迭代次数
	iter_num = iter_num + 1;
	S = loadbus;
	S_to = S;
	S_fr = S_to;
	% 前推过程
	    for k = NUMLINE:-1:2                     % 从最后一条支路,依次-1, 至第一条支路
		i = firstbranchbus(k);          % 每个支路的首节点
		S_fr(k) = S_to(k) + Zbranch(k) * abs(S_to(k)/U(k))^2;  
		S_to(i) = S_to(i) + S_fr(k);
        end

	% 回推过程

	    for k = 2:NUMLINE % 从第一条支路,依次+1, 至最后一条支路
		i = firstbranchbus(k);
		U(k) = U(i) - Zbranch(k) * conj(S_fr(k)/U(i));
        end

	% 检查收敛

	dU = abs(U - Uold); % 偏差量
	Uold = U;


    %PV
	if ~isempty(PVNUM) % 存在PV节点 
	% 计算PV节点的Q校正   假设P=0
        Ug = Upv(length(PVNUM));                            % 读所有PV电源控制的电压值
		dE = ( Ug./abs(U(PVNUM))-1).*real(U(PVNUM) );      %电压差值
		dD = -Bpv * dE;                                      %补偿电流虚部
		dC = dD .* imag(U(PVNUM))./real(U(PVNUM));           %补偿电流实部
	    V_compen = GenVbuchang(dC+1j*dD, PVNUM, NUMBUS, NUMLINE, firstbranchbus, Zbranch); % 生成节点电压补偿量
		U = U + V_compen;                                      % 节点电压迭代
		DQ = -dD .* abs(U(PVNUM)).^2 ./ real(U(PVNUM));        % PV节点 无功输出量
		Qpv = Qpv + DQ;                                        % 更新 PV节点 无功功率输出
		loadbus(PVNUM) = loadbus(PVNUM) - 1j*DQ;               % 将PV节点输出的无功 代入节点功率中用于计算
	end
	

	if  ~isempty(PINUM) % 存在PI节点
		Q_PI = sqrt( Ipi_G^2 * abs(real(U(PINUM))^2 + imag(U(PINUM))^2) - Ppi_G^2); %该这么求吗？
		if abs(Q_PI) > Ppi_G % 判断无功输出是否越界，一般无功输出不能超过有功输出范围，可视为有功输出最大值的正负范围内，即 QG = [-PG,+PG]
			Q_PI = Ppi_G;
		end
		loadbus(PINUM) = loadbus(PINUM) - 1j*Q_PI;  % 将PI节点计算得出的无功 代入节点功率中用于计算
		
    end

	if max(dU) <= ep    % 若达到收敛精度
		Uppf=Uppf+(1/NMC).*U;                     %求MC的期望值
        Uprintbus(mci,:)=U(:,1);  %求所有节点每次的模拟值。
        for t=1:NUMBUS
             if ( (abs(U(t))<0.93) || (abs(U(t))>1.08)   )  %电压越限
                    Ucrossline(t,1)= Ucrossline(t,1)+1;
             end
        end 

        %Sslack = S_to(1);
        S_fr = S_fr(2:end); % i->j功率
        S_to = S_to(2:end); % j->i功率
        Sloss = S_fr - S_to;
        Ploss = SB*real(Sloss); % kW 有功网损
        Qloss = SB*imag(Sloss); % kVar 无功网损
        Sloss = Ploss+1j*Qloss;
     
        plossppf(1,mci)=sum(Ploss);
        break
    end

       
    end

end



exetime= etime(clock, time0);  %程序运行的时间
 

% 保存发电机无功输出 和 支路功率流
if ~isempty(PVNUM)
	Qpv_G = Qpv; % p.u.
end





%% 进行电压概率分布图绘制

fprintf('程序计算时间为%.2f seconds\n', exetime);

%图-电压频率直方图
k=100;
figure(4);
histogram(abs(Uprintbus(:,printnumbus)),100,'Normalization','probability');
hold on
histogram(abs(RANDOMMC.UMC(:,printnumbus)),100,'Normalization','probability');
title('电压频率分布直方图');xlabel('Vppf');ylabel('概率');  %建立这一节点电压概率模型。
legend('LHS采样','MC随机采样')
%图-电压经验累计分布图CDF
%cdfplot(abs(Uprintbus(:,printnumbus)));
%title('电压经验分布曲线');xlabel('Vppf');ylabel('F(x)');

%% 画网损概率图

figure(5)
h1 = histogram(plossppf,200,'Normalization','probability');

title('有功网损分布图');xlabel('Ploss(kW)');ylabel('概率');
%% 画各节点电压数据图

 figure(2)
 plot(angle(U),'-k<')
 title('各节点相位')
 xlabel('节点编号')
 ylabel('节点相位（°）')
 grid on
 hold on

 figure(1)
 plot(abs(U),'-k<')
 hold on
 plot(abs(Uppf),'-r<')
 hold on
 plot(abs(RANDOMMC.UPPF),'-g<')
 hold on
 title('各节点电压幅值')
 xlabel('节点编号')
 ylabel('节点电压')
 grid on
 hold on
 legend('非概率潮流','LHS采样概率潮流','MC随机抽样概率潮流')

figure(2)
 plot(angle(Uppf),'-r<')
 title('各节点相位')
 xlabel('节点编号')
 ylabel('节点相位（°）')
 grid on
 hold on
 legend('非概率潮流','MC概率潮流')

 figure(3)
 plot(Ucrossline/NMC,'-k<')
 title('各节点电压越限概率')
 xlabel('节点编号')
 ylabel('概率(%)')
 grid on
 hold on

 bianyi=RANDOMMC.BIANYI;

 figure(6)
 plot(RANDOMMC.BIANYI,'-g<')
 hold on
 plot(LHSMCMC.BIANYI,'-k<')
 title('各节点电压计算精度')
 xlabel('节点编号')
 ylabel('变异系数')
 grid on
 hold on
 legend('随机模拟MC','LHSMC')