function [ ] = FOBD_implicit(CFL0,deltaX0,ui0,tf0,a0,N0,NP0)
% Function to Solve Wave equation using First order backward difference,
% implicit integration. Inputs to the function is CFL0 number, Grid Size,
% Initial condition, Final time, speed of wave, Number of Mesh points and
% Plotting location in terms of number of time steps

deltaT=CFL0*deltaX0/a0;
t0=0;
tplot=0;
usol=ui0;
usoln=usol;
ufinal=zeros(N0,1);
ufinal(2:N0-1)=usoln;
ufinal(1,1)=0;
ufinal(N0,1)=0;

a=-CFL0;
b=1+CFL0;

fig=figure('Name',['First-order backward difference, implicit integration for CFL=',num2str(CFL0)]);
title(['First-order backward difference implicit integration for CFL=',num2str(CFL0)]);
xlabel('x')
ylabel('u')
hold on;
umax=0;
umin=0;

while(t0-tf0<1e-6)
%Plotting
if(abs(t0-tplot)<1e-5)
tplot=t0+NP0*deltaT;
x=linspace(0,1,N0);
plot_customized(x,ufinal,t0,'b');
umax=max(umax,max(ufinal));
umin=min(umin,min(ufinal));
end

%Finding the Solution
usoln(1)=usol(1)/b;
for i=2:N0-2
usoln(i)=(usol(i)-a*usoln(i-1))/b;
end

usol=usoln;
ufinal(2:N0-1)=usoln;
t0=t0+deltaT;
end
axis([0 1 umin-0.1 umax+0.15])
str=['6FOBD_implicit_integration_for_CFL_',num2str(CFL0),'.png']
saveas(fig,str);
end

