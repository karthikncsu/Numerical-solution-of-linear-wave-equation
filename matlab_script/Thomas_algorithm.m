function [ usol ] = Thomas_algorithm(K,f)
% Function uses thomas algorithm to Solve Tridiagonal System of Matrices
[N1 N2]=size(K);
usol=zeros(N1,1);
% Eliminating right Diagnol Elements    
for i=2:N1                     
K(i,i)=K(i,i)-K(i-1,i)/K(i-1,i-1)*K(i,i-1);
f(i)=f(i)-f(i-1)/K(i-1,i-1)*K(i,i-1);
K(i,i-1)=0;
end

usol(N1)=f(N1)/K(N1,N1);   % Finding u-values
for i=N1-1:-1:1
usol(i)=(f(i)-usol(i+1)*K(i,i+1))/K(i,i);
end


end

