function [ ] = FourthOBD_implicit(CFL0,deltaX0,ui0,tf0,a0,N0,L0,NP0)
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

a=CFL0/3;
b=-3*CFL0/2;
c=3*CFL0;
d=1-11*CFL0/6;


% %Matrix Assebling
% K=zeros(N0-2,N0-2);
% for i=4:1:N0-2
%     K(i,i-3)=-CFL0/3;
%     K(i,i-2)=3*CFL0/2;
%     K(i,i-1)=-3*CFL0;
%     K(i,i)=1+11*CFL0/6;
%  end
% K(1,1)=1+11*CFL0/6+3*CFL0/2;
% K(1,2)=-CFL0/3;
% K(2,1)=-CFL0/3-3*CFL0;
% K(2,2)=1+11*CFL0/6;
% K(3,1)=3*CFL0/2;
% K(3,2)=-3*CFL0;
% K(3,3)=1+11*CFL0/6;

%Creating Figure
fig=figure('Name',['Fourth order backward difference , imlicit integration for CFL=',num2str(CFL0)]);
title(['Fourth order backward difference, imlicit integration for CFL=',num2str(CFL0)]);
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
%Finding the Solution
for i=4:N0-2
usoln(i)=a*usol(i-3)+b*usol(i-2)+c*usol(i-1)+d*usol(i);
end
usoln(1)=(b+d)*usol(1)+a*usol(2);
usoln(2)=(a+c)*usol(1)+d*usol(2);

%usoln=K\(usol);
% f=K1*usol;
% usoln=Thomas_algorithm(K,f);
usol=usoln;
ufinal(2:N0-1)=usoln;
t0=t0+deltaT;
end
% plot_customized(x,ufinal,t0,'r'); %Plotting final Solution
axis([0 1 umin-0.1 umax+0.15])
str=['5FourthOrder_backward_difference_imlicit_for_CFL_',num2str(CFL0),'.png']
saveas(fig,str);

end