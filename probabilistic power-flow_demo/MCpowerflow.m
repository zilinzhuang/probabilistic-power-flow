clc
clear
close all
%% 含DG、随机模拟抽样的概率潮流计算
%% 假设条件

MODEL = data33origin; % 导入33节点数据。负荷、支路阻抗数据均为标幺值

SB = 10;        % 10MVA 基准功率
SB = SB*1000;   % 基准功率转换为kVA
UB = 12.66;     % 基准电压 kV
IB = SB/(sqrt(3)*UB); % 基准电流 A


U_ref = 1.0; % 首节点/平衡节点
NUMBUS = 33; NUMLINE = 32; % 总节点数, 总线路数
firstbranchbus = MODEL.branch(:,2); % 首节点列向量
Zbranch = MODEL.branch(:,4)+1i*MODEL.branch(:,5); % 支路阻抗矩阵p.u.
Ybranch = 1./Zbranch;                         % 支路的电纳矩阵p.u.
loadbus = MODEL.bus(:,3)+1i*MODEL.bus(:,4);           % 负荷阻抗矩阵 p.u.
Sd_orginal = loadbus;                             % 原负荷数据
ep = 1e-5; 
iter_max = 100;                 % 收敛精度及最大迭代次数
U = U_ref * ones(NUMBUS,1);
Uold = U;
iter_num = 0;



firstbranchbus = [0; firstbranchbus]; % 添加一条虚拟支路至所有支路的顶层, 这样支路数就等于节点数，便用于后面迭代计算
Zbranch = [0; Zbranch];    %假设与首节点连接的支路阻抗为0，并把支路数加1.
NUMLINE = NUMLINE + 1;             
Uppf=zeros(NUMBUS,1);                     %概率潮流电压
Ucrossline=zeros(NUMBUS,1);                % 电压越限数
Slossppf=zeros(NUMLINE-1,1);              %概率潮流下的线路损耗
%% Monte-Carlo概率模型建立

NMC=2000;                          %  monte-carlo模拟数

%% 随机抽样MC
wind1=wblrnd(7,2,1,NMC);                %  设定WEIBULL参数！scale parameter A and shape parameter B.后两项为矩阵形式
sun1=betarnd(4,0.1,1,NMC);             %  设定BETA参数！
pw1=wind2power(wind1,5,1)*1000;       %  风能有功功率
ps1=solar2power(sun1,700)*1000;       %  太阳能功率




fprintf('请输入需要哪个节点的电压概率图像？(范围1-33)\n');
printnumbus=input(['节点号=']);
Uprintbus=zeros(NMC,NUMBUS);          %建立所有模拟的数据库矩阵

plossppf=zeros(1,NMC);                %有功网损
time0 = clock;  % 计时开始

%% 进行广度优先搜索编号

stemp=[1];      %这一层搜索储存的节点号
sstemp=[];      %循环中储存的节点
n=1;
m=1;
TempBranch=MODEL.branch;  %目前可利用的支路（去除已经存的支路）

while ~isempty(TempBranch)
    for k=1:length(stemp)
        for i=1:size(TempBranch,1)   
            if TempBranch(i,2)==stemp(k)
            s1(n,:)=TempBranch(i,:);            %储存上一层节点连通的支路，并记录连接的节点号
            n=n+1;
            sstemp=[sstemp,TempBranch(i,3)];
            else
            s2(m,:)=TempBranch(i,:);            %若有不满足的支路，那么继续存在矩阵内。
            m=m+1;
            end

        end
    TempBranch=s2;           %将s2赋值给TempBranch 重复上述判断，直到TempBranch为空集为止
    s2=[];
    m=1;

    end
stemp=sstemp;                %这一次的临时变量，成为下一次的初始值（首节点）
sstemp=[];
end
s3=s1;
for(i=1:size(s1,1))
    s1(i,:)=s3(size(s1,1)-i+1,:);       %逆序排列
end



for mci=1:NMC 

%初始化

loadbus = Sd_orginal;    % 负荷阻抗矩阵 p.u.
U = U_ref * ones(NUMBUS,1);
Uold = U;
iter_num = 0; 

%% 对DG分类、定义

% PQ型 DG        双馈风机、同步全功率风机
PQNUM = []; % PQ型DG位置          
Ppq_G = [pw1(1,mci)/SB , pw1(1,mci)/SB]; % kW >>> p.u.
Qpq_G = Ppq_G.*tan(acos(0.98)); % kVar >>> p.u  %假定功率因数恒定

% PI型 DG
PINUM = [31];    % PI型DG位置、光伏（电流型逆变器）        
Ppi_G = [ps1(1,mci)/SB]; % kW >>> p.u.
Ipi_G = [-10/IB]; % A >>> p.u
Q_PI=zeros(1,length(PINUM));
Q_OLDPI=zeros(1,length(PINUM));


% PQ(V)型 DG
PQVNUM = [18];    % PQV型DG位置  定速型异步风机（转差率不变，导致电压不变、无功吸收不变）
Ppqv_G = [30/SB]; % kW >>> p.u.
s_pqv=0.3;       % 转差率
Qpqv_G = -Ppqv_G.*(100^2+20*(250+20)*(s_pqv^2)) / (100*250*s_pqv);      % 求解建立磁场的无功吸收量

% PV型 DG         燃料电池、光伏（电压型逆变器）
PVNUM = [13]; % PV型DG位置
Ppv_G = [100/SB]; % kW >>> p.u
Qpv_G = [200/SB]; % PV型DG的无功出力调节范围，这里仍按照恒功率模型计算。当无功大于此值，转为PQ节点。kVar>>> p.u
Upv = [0.98]; % PV型DG 定电压
npv = length(PVNUM);
Qpv = zeros(npv,1);           % 建立PV无功矩阵


%% 负荷加入概率模型的变化，接PQ、PQV、PI、PV节点的更新。
Loadkplf=normrnd(1,0.1,33,1); %生成33*1的正态分布的随机数矩阵。
loadbus(:,1)=Sd_orginal(:,1).*Loadkplf(:,1);

for t=1:length(PQNUM)
loadbus(PQNUM(t)) = loadbus(PQNUM(t)) - complex(Ppq_G(t),Qpq_G(t));  % 对PQ DG 把其作为负的负荷加到节点负荷处
end


for t=1:length(PINUM)
loadbus(PINUM(t)) = loadbus(PINUM(t)) - Ppi_G(t);               % 对PI DG 先将其恒定的有功输出作为负的负荷添加到节点负荷处
end


for(t=1:length(PVNUM))
loadbus(PVNUM(t)) = loadbus(PVNUM(t)) - Ppv_G(t);               % 对PV DG 先将其恒定的有功输出作为负的负荷添加到节点负荷处
end

for(t=1:length(PQVNUM))
loadbus(PQVNUM(t)) = loadbus(PQVNUM(t)) - complex(Ppqv_G(t),Qpqv_G(t));               % 对PQV DG 把其作为负的负荷加到节点负荷处
end



%% 前推回代潮流
            

while iter_num < iter_max             %迭代次数
	iter_num = iter_num + 1;
	S = loadbus;
	S_to = S;
	S_fr = S_to;
	% 前推过程
	    for k = 1:1:NUMLINE-1                     % 从排序好的第一条支路，开始+1
            stemp=s1(k,:);
		head = stemp(2);          % 每个支路的首节点
        tail = stemp(3);           % 尾节点
		S_fr(tail) = S_to(tail) + (stemp(4)+ 1i*stemp(5))* abs(S_to(tail)/U(tail))^2;  
		S_to(head) = S_to(head) + S_fr(tail);
        end

	% 回推过程

	    for k = NUMLINE-1:-1:1 % 从排序好的最后支路,依次-1, 至第一条支路
		stemp=s1(k,:);
        head = stemp(2);          % 每个支路的首节点
        tail = stemp(3);          % 每个支路的尾节点
		U(tail) = U(head) -(stemp(4)+ 1i*stemp(5)) * conj(S_fr(tail)/U(head));
        end

	% 计算偏差量

	dU = abs(U - Uold); % 偏差量
	Uold = U;

%% PI节点处理        

	if  ~isempty(PINUM)                % 存在PI节点
        for s=1:length(PINUM)
		Q_PI(1,s) = sqrt (Ipi_G(1,s)^2 * abs(U(PINUM(1,s))).^2  - Ppi_G(1,s)^2); %根据给定的P、I参数，算出该节点的无功功率，转换为PQ节点。
		    if abs(Q_PI(1,s)) > Ppi_G(1,s) % 判断无功输出是否越界.这里假定无功不能超过有功输出范围，最大无功输出为Ppi_G
			Q_PI(1,s) = Ppi_G(1,s);
		    end
		loadbus(PINUM(1,s)) = loadbus(PINUM(1,s)) - 1j*(Q_PI(1,s)-Q_OLDPI(1,s));  % 将PI节点计算得出的无功 代入节点功率中用于计算
        end
        Q_OLDPI=Q_PI;
    end


      %% PV节点处理
    % 生成PV节点的阻抗矩阵 为了求解每增加单位电流下，该节点电压增量是多少。
    Zpv = GenZpv(PVNUM,NUMBUS,NUMLINE,firstbranchbus,Zbranch); % 子函数 
    if (size(Zpv,1) > 0)
	Xpv = (imag(Zpv));              % 电抗，以便后续利用
    end

	if ~isempty(PVNUM) % 存在PV节点 
	% 计算PV节点的Q校正   假设ΔP=0
        for s=1:npv
        Ug = Upv(1,s);                            % 读所有PV电源控制节点的电压值（设定值）
		dE = ( Ug./abs(U(PVNUM(1,s)))-1).*real(U(PVNUM(1,s)) );        %电压差值
		dD = -dE./Xpv(s);                                            %补偿电流虚部
		dC = dD .* imag(U(PVNUM(1,s)))./real(U(PVNUM(1,s)));           %补偿电流实部
	    V_compen = GenVbuchang(dC+1j*dD, PVNUM(1,s), NUMBUS, NUMLINE, firstbranchbus, Zbranch); % 通过补偿电流，前推回代，生成每个节点在这个电流作用下的电压补偿量
		U = U + V_compen;                                      % 节点电压更新
		DQ = -dD .* abs(U(PVNUM(1,s))).^2 ./ real(U(PVNUM(1,s)));        % PV节点 无功输出量
		Qpv(s,1) = Qpv(s,1) + DQ;
        
        if Qpv_G(1,s)> Qpv(s,1)                                       % 判断无功功率是否越限
         loadbus(PVNUM(1,s))=Sd_orginal(PVNUM(1,s))-Qpv_G(1,s);       %若越限，PV节点转为PQ节点
        end
            
            % 更新 PV节点 无功功率输出
		loadbus(PVNUM(1,s)) = loadbus(PVNUM(1,s)) - 1j*DQ;               % 将PV节点输出的无功 代入节点功率中用于计算
	    end
	
    end
	

    %% 判断收敛性

	if max(dU) <= ep    % 若达到收敛精度
	
        Uppf=Uppf+(1/NMC).*U;                     %求MC的期望值
        Uprintbus(mci,:)=U(:,1);  %求所有节点每次的模拟值。
        for t=1:NUMBUS
             if ( (abs(U(t))<0.93) || (abs(U(t))>1.07)   )  %电压越限
                    Ucrossline(t,1)= Ucrossline(t,1)+1;
             end
        end 
      
        S_fr = S_fr(2:end); % i>>j功率
        S_to = S_to(2:end); % j>>i功率
        Sloss = S_fr - S_to;
        Ploss = SB*real(Sloss); % kW 有功网损
        Qloss = SB*imag(Sloss); % kVar 无功网损
        Sloss = Ploss+1j*Qloss;
        Slossppf(:,1)=Sloss(:,1)+Slossppf(:,1);
        plossppf(1,mci)=sum(Ploss);
        break
    end

       
    end

end





exetime= etime(clock, time0);  %程序运行的时间
 







%% 进行电压概率分布图绘制

fprintf('程序计算时间为%.2f seconds\n', exetime);

%图-电压频率直方图
k=100;
figure(4);
histogram(abs(Uprintbus(:,printnumbus)),100,'Normalization','probability');
hold on
title('电压频率分布直方图');xlabel('Vppf');ylabel('概率');  %建立这一节点电压概率模型。
legend('MC随机采样')
%图-电压经验累计分布图CDF
%cdfplot(abs(Uprintbus(:,printnumbus)));
%title('电压经验分布曲线');xlabel('Vppf');ylabel('F(x)');

%% 画网损概率图

figure(5)
h1 = histogram(plossppf,200,'Normalization','probability');

title('有功网损分布图');xlabel('Ploss(kW)');ylabel('概率');

figure(7)
 plot(imag(Slossppf/NMC),'-g<')
 hold on
 plot(real(Slossppf/NMC),'-k<')
 title('无功线损和有功线损')
 xlabel('节点编号')
 ylabel('功率')
 grid on
 hold on
 legend('无功线损','有功线损')

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

 title('各节点电压幅值')
 xlabel('节点编号')
 ylabel('节点电压')
 grid on
 hold on
 legend('非概率潮流','MC随机抽样概率潮流')

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

 

 figure(6)
 plot(std(abs(Uprintbus))./mean(abs(Uprintbus)) ,'-g<')
 hold on
 title('各节点电压计算精度')
 xlabel('节点编号')
 ylabel('变异系数')
 grid on
 hold on
 legend('随机模拟MC')