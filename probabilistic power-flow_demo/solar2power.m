function ps=solar2power(solar,n)
%%  根据光强的抽样样本换算出光伏电站的有功出力样本,solar为光强向量，n为太阳能光伏电池板的个数
    A=4;    %每块板的面积
    g=0.25;  %效率
    solar(find(solar>1))=1;
     ps=solar*A*g*n/1000;
 
     
     return;
 
     