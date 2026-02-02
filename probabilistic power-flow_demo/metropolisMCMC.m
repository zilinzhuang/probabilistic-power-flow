
d =zeros(1,400000);
x=2;
for i=1:400000
    y = unifrnd(x-1,x+1);
    if y<0
        y=x;
    end
    h=min(1, wblpdf(7,2,y)/ wblpdf(7,2,x));
    U=rand;           %产生均值为0，方差 σ^2 = 1，标准差σ = 1的正态分布的随机数或矩阵的函数。
    if U<h
        x=y;
    end
    d(i)=x;
end
a=0:0.08:20;
hist(d(200000:400000),a)
