function pw=wind2power(vw,n,sty)
%%  根据风速的抽样样本换算出风机的有功出力样本
%% vw为风速的抽样样本，n为风电场装设风电机组的台数，sty为风电机组的类型
 %选择：
 %           1   异步风电机组
 %           2   双馈风电机组 DFIG
 %           3   PMSG
 % 三种机型对应的风速功率曲线不同（怎么详细确定它的功率曲线？）
 
 if sty==1
     x=[0:50];
     y=[0 0 0 0 ,0.5/16*[4:20]-0.125,0.5*ones(1,10),zeros(1,20)];%why?
     pw=interp1(x,y,vw);           %线性插值，用spline会不会比linear更好？
     pw=pw*n;
 else
     if sty==2
     x=[0:50];
     y=[0 0 0 0 0,0.05*[5:15]-0.25,0.5*ones(1,10),zeros(1,25)];
     pw=interp1(x,y,vw);
     pw=pw*n;
     end
     if sty==3
     x=[0:50];
     y=[0 0 0 0 0,0.05*[5:15]-0.25,0.5*ones(1,5),(-0.5/15)*[21:36]+(18/15),zeros(1,14)];
     pw=interp1(x,y,vw);
     pw=pw*n;
     end
     return;
 end