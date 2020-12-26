% Program to Solve Linear Wave Equation
% ut+a*ux=0
clc
clear all
close all
a0=1; %Speed of Wave
L0=1; % Length of the Interval
deltaX=0.005; %Grid Size
CFL=[0.4; 1; 1.3]; % Required CFL Numbers to run
N=L0/deltaX+1; %Number of grid points in the interval
xf=0.75; %required position of the discontinuity
NP=40; %Number of timesteps to plot the solution

%Initial Conditions
ui=zeros(N-2,1);

for i=1:1:N-2
 x=(i-1)*deltaX;
% ui(i)=sin(2*pi*10*x/L0);
    if(x>0.2&&x<=0.3)
     ui(i)=1;   
    end
end

tf=(xf-0.3)/a0;
 
%Calling First Order Backward Difference Method
  FOBD_explicit(CFL(1),deltaX,ui,0.45,a0,N,round(tf/5/(CFL(1)*deltaX/a0)));
  FOBD_explicit(CFL(2),deltaX,ui,0.45,a0,N,round(tf/2/(CFL(2)*deltaX/a0)));
  FOBD_explicit(CFL(3),deltaX,ui,0.2,a0,N,round(0.1/1/(CFL(2)*deltaX/a0)));

%Calling Second Order Backward Difference Method
SOBD_explicit(CFL(1),deltaX,ui,0.1,a0,N,round(tf/10/(CFL(1)*deltaX/a0)));  
SOBD_explicit(CFL(2),deltaX,ui,0.08,a0,N,round(tf/10/(CFL(2)*deltaX/a0)));
SOBD_explicit(CFL(3),deltaX,ui,0.08,a0,N,round(tf/10/(CFL(3)*deltaX/a0)));

% Calling Crank-Nicolson Implicit Method
CN_implicit(CFL(1),deltaX,ui,0.2,a0,N,round(tf/5/(CFL(1)*deltaX/a0)));  
CN_implicit(CFL(2),deltaX,ui,0.2,a0,N,round(tf/5/(CFL(2)*deltaX/a0))); 
CN_implicit(CFL(3),deltaX,ui,0.2,a0,N,round(tf/5/(CFL(3)*deltaX/a0)));
 
%Calling Lax-Wendroff Scheme
LW_explicit(CFL(1),deltaX,ui,0.45,a0,N,round(0.45/5/(CFL(1)*deltaX/a0)));  
LW_explicit(CFL(2),deltaX,ui,0.45,a0,N,round(0.45/2/(CFL(2)*deltaX/a0))); 
LW_explicit(CFL(3),deltaX,ui,0.1,a0,N,round(0.1/1/(CFL(3)*deltaX/a0)));

%Calling Lax Scheme
LAX_explicit(CFL(1),deltaX,ui,0.45,a0,N,round(tf/5/(CFL(1)*deltaX/a0)));  
LAX_explicit(CFL(2),deltaX,ui,0.45,a0,N,round(tf/2/(CFL(2)*deltaX/a0)));  
LAX_explicit(CFL(3),deltaX,ui,0.2,a0,N,round(tf/5/(CFL(3)*deltaX/a0)));  

% Calling First order implicit scheme
FOBD_implicit(CFL(1),deltaX,ui,0.45,a0,N,round(0.45/5/(CFL(1)*deltaX/a0)));  
FOBD_implicit(CFL(2),deltaX,ui,0.45,a0,N,round(0.45/5/(CFL(2)*deltaX/a0)));  
FOBD_implicit(CFL(3),deltaX,ui,0.45,a0,N,round(0.45/5/(CFL(3)*deltaX/a0)));  



