function [ ] = LW_explicit(CFL0,deltaX0,ui0,tf0,a0,N0,NP0)
% Function to Solve Wave equation using Lax Wendroff explicit scheme
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

a=0.5*CFL0*(1+CFL0);
b=1-CFL0*CFL0;
c=0.5*CFL0*(CFL0-1);

fig=figure('Name',['Lax-Wendroff Explicit Scheme for CFL=',num2str(CFL0)]);
title(['Lax-Wendroff Explicit Scheme for CFL=',num2str(CFL0)]);
xlabel('x')
ylabel('u')
hold on;
umax=0;
umin=0;
while(t0-tf0<1e-6)
%Plotting after NP0 number of time steps
if(abs(t0-tplot)<1e-5)
tplot=t0+NP0*deltaT;
x=linspace(0,1,N0);
plot_customized(x,ufinal,t0,'b');
umax=max(umax,max(ufinal));
umin=min(umin,min(ufinal));
end

%Finding the Solution
for i=2:N0-3
usoln(i)=a*usol(i-1)+b*usol(i)+c*usol(i+1);
end
usoln(1)=b*usol(1)+c*usol(2);
usoln(N0-2)=a*usol(N0-3)+b*usol(N0-2);

usol=usoln;
ufinal(2:N0-1)=usoln;
t0=t0+deltaT;
end

%plot_customized(x,ufinal,t0,'r'); %Plotting final Solution
axis([0 1 umin-0.1 umax+0.15])
str=['4Laxwendroff_imlicit_integration_for_CFL_',num2str(CFL0),'.png'];
saveas(fig,str);
end