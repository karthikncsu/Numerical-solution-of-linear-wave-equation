function [ ] = CN_implicit(CFL0,deltaX0,ui0,tf0,a0,N0,NP0)
% Function to Solve Wave equation using Crank-Nicolson implicit scheme
% Inputs to the function is CFL0 number, Grid Size, Initial condition,
% Final time, speed of wave, Number of Mesh points Plotting location in
% terms of number of time steps

deltaT=CFL0*deltaX0/a0;
t0=0;
tplot=0;
usol=ui0;
usoln=usol;
ufinal=zeros(N0,1);
ufinal(2:N0-1)=usoln;
ufinal(1,1)=0;
ufinal(N0,1)=0;

%Matrix Assebling
K=zeros(N0-2,N0-2);
K1=zeros(N0-2,N0-2);
for i=2:1:N0-3
    K(i,i-1)=-CFL0/4;
    K(i,i)=1;
    K(i,i+1)=CFL0/4;
    K1(i,i-1)=CFL0/4;
    K1(i,i)=1;
    K1(i,i+1)=-CFL0/4;
end
K(1,1)=1;
K(1,2)=CFL0/4;
K(N0-2,N0-3)=-CFL0/4;
K(N0-2,N0-2)=1;
K1(1,1)=1;
K1(1,2)=-CFL0/4;
K1(N0-2,N0-3)=CFL0/4;
K(N0-2,N0-2)=1;

%Creating Figure
fig=figure('Name',['Crank-Nicolson, imlicit integration for CFL=',num2str(CFL0)]);
title(['Crank-Nicolson, imlicit integration for CFL=',num2str(CFL0)]);
xlabel('x')
ylabel('u')
hold on;
umax=0;
umin=0;

while(t0-tf0<1e-6)
%Plotting after N0 number of time steps
x=linspace(0,1,N0);
if(abs(t0-tplot)<1e-5)
tplot=t0+NP0*deltaT;
plot_customized(x,ufinal,t0,'b');
umax=max(umax,max(ufinal));
umin=min(umin,min(ufinal));
end

f=K1*usol;
usoln=Thomas_algorithm(K,f);
usol=usoln;
ufinal(2:N0-1)=usoln;
t0=t0+deltaT;
end
% plot_customized(x,ufinal,t0,'r'); %Plotting final Solution
axis([0 1 umin-0.1 umax+0.15])
str=['3Crank_Nicolson_imlicit_integration_for_CFL_',num2str(CFL0),'.png']
saveas(fig,str);

end