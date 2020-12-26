function [ ] = CTCD_explicit(CFL0,deltaX0,ui0,tf0,a0,N0,L0,NP0)
% Function to Solve Wave equation using First order backward difference
% Inputs to the function is CFL0 number, Grid Size, left and right boundary
% conditions, Inital Conditions, Final time, Speed of the wave
deltaT=CFL0*deltaX0/a0;
t0=0;
tplot=0;
usol=ui0;
usoln=usol;
usoln0=ui0;
% usoln0=zeros(N0-2,1);
ufinal=zeros(N0,1);
ufinal(2:N0-1)=usoln;
ufinal(1,1)=0;
ufinal(N0,1)=0;

%Matrix Assebling
K=zeros(N0-2,N0-2);
for i=2:1:N0-3
    K(i,i-1)=CFL0;
    K(i,i+1)=-CFL0;
end
K(1,2)=-CFL0;
K(N0-2,N0-3)=CFL0;

size(K)
size(usoln0)
size(usol)
size(usoln)

fig=figure('Name',['Central Time Central difference for CFL=',num2str(CFL0)]);
title(['Central Time Central Difference for CFL=',num2str(CFL0)]);
xlabel('x')
ylabel('u')
hold on;
umax=0;
umin=0;

while(t0-tf0<1e-6)
%Plotting
if(abs(t0-tplot)<1e-5)
tplot=t0+NP0*deltaT;
x=linspace(0,L0,N0);
plot_customized(x,ufinal,t0,'b');
umax=max(umax,max(ufinal));
umin=min(umin,min(ufinal));
end

%Finding Solution
usoln=K*usol+usoln0;
usoln0=usol;
usol=usoln;
ufinal(2:N0-1)=usoln;
t0=t0+deltaT;
end
axis([0 1 umin-0.1 umax+0.15])
str=['CTCD_explicit_integration_for_CFL_',num2str(CFL0),'.png']
saveas(fig,str);
end