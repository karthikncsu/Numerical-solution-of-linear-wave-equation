function [ ] = Fourth_implicit(CFL0,deltaX0,ui0,tf0,a0,N0,L0,NP0)
% Function to Solve Wave equation using Fourth order central difference implicit scheme. Inputs
%to the function is CFL0 number, Grid Size, left and right boundaryconditions, Inital Conditions,
%Final time, Speed of the wave, Plotting location in terms on number of time steps
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
for i=3:1:N0-4
    K(i,i-2)=CFL0/12;
    K(i,i-1)=-2*CFL0/3;
    K(i,i)=1;
    K(i,i+1)=2*CFL0/3;
    K(i,i+2)=-CFL0/12;
 end
K(1,1)=1+CFL0/12;
K(1,2)=2*CFL0/3;
K(1,3)=-CFL0/12;
K(2,1)=-2*CFL0/3;
K(2,2)=1;
K(2,3)=2*CFL0/3;
K(2,4)=-CFL0/12;

K(N0-3,N0-5)=CFL0/12;
K(N0-3,N0-4)=-2*CFL0/3;
K(N0-3,N0-3)=1;
K(N0-3,N0-2)=2*CFL0/3;
K(N0-2,N0-4)=CFL0/12;
K(N0-2,N0-3)=-2*CFL0/3;
K(N0-2,N0-2)=1-CFL0/12;

%Creating Figure
fig=figure('Name',['Fourth order central difference , imlicit integration for CFL=',num2str(CFL0)]);
title(['Fourth order central difference, imlicit integration for CFL=',num2str(CFL0)]);
xlabel('x')
ylabel('u')
hold on;
umax=0;
umin=0;

while(t0-tf0<1e-6)
%Plotting after N0 number of time steps
x=linspace(0,L0,N0);
if(abs(t0-tplot)<1e-5)
tplot=t0+NP0*deltaT;
plot_customized(x,ufinal,t0,'b');
umax=max(umax,max(ufinal));
umin=min(umin,min(ufinal));
end

%Finding Solution
usoln=inv(K)*(usol);
% f=K1*usol;
% usoln=Thomas_algorithm(K,f);
usol=usoln;
ufinal(2:N0-1)=usoln;
t0=t0+deltaT;
end
% plot_customized(x,ufinal,t0,'r'); %Plotting final Solution
axis([0 1 umin-0.1 umax+0.15])
str=['5FourthOrder_central_difference_imlicit_for_CFL_',num2str(CFL0),'.png']
saveas(fig,str);

end