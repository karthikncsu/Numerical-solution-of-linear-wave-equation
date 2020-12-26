function [ ] = FOBD_explicit(CFL0,deltaX0,ui0,tf0,a0,N0,NP0)
% Function to Solve Wave equation using First order backward difference,
% explicit integration Inputs to the function is CFL0 number, Grid Size,
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

b=1-CFL0;
a=CFL0;

fig=figure('Name',['First-order backward difference, explicit integration for CFL=',num2str(CFL0)]);
title(['First-order backward difference explicit integration for  CFL=',num2str(CFL0)]);
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
for i=2:N0-2
usoln(i)=a*usol(i-1)+b*usol(i);
end
usoln(1)=b*usol(1);

usol=usoln;
ufinal(2:N0-1)=usoln;
t0=t0+deltaT;
end
axis([0 1 umin-0.1 umax+0.15])
str=['1FOBD_explicit_integration_for_CFL_',num2str(CFL0),'.png'];
saveas(fig,str);
end

