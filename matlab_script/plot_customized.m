function [  ] = plot_customized(x,ufinal,t0,col)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

plot(x,ufinal,'color',col,'LineWidth',2);
hold on 
grid on;
[max_val I]=max(ufinal);

%axis([0 1 -0.25 1.25])

str = ['t=',num2str(t0)];
text(x(I),max_val+0.05,str);
hold on


end

