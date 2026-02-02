
clc
clear all
close all
a = randn(10000,1);
h1 = histogram(a,21,'Normalization','pdf');
hold on
x = -4:0.1:4;
mu = 0;
sigma = 1;
y = exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(x,y,'LineWidth',1.5)
figure(2)
h2 = histogram(a,21);
L = length(h2.BinEdges);
for i = 1:L-1
    x2(i) = (h2.BinEdges(i) + h2.BinEdges(i+1))/2;
end
y2 = h2.Values/10000/h2.BinWidth;
figure(1)
hold on
plot(x2,y2,'o','linewidth',1.5,'MarkerSize',8)
xlim([-4,4])

