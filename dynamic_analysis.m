%%% the MATLAB file is in the attachment 
clear all
close all 
clc
%%%% this code is to find the fundamentale natural frequency for a mast
%%%% carrying a horizontal axis wind turbine 
%%%the 1 st section is to calculate the area and define the cantiliver beam
%%%
rho=8050;                                %the mass denisty of the mast 
L=60;                                    % the Height of the mast 
D_og=4;   r_og=D_og/2;                   %the outer diameter of the lower part of the mast
D_ou=2;   r_ou=D_ou/2;                   %the outer diameter of the lower part of the mast         
t=.0508;                                 % the thickness of the mast 
e=10;                                    % number of elements in the global system
n1=(2*e)+2;                              % the number of nodes in the global system 
NI=e+1;
Lx=linspace(0,L,NI);
Li=Lx(2)-Lx(1);                          %the element lengtth
ri=r_og-t;  
A=pi*(r_og.^2-ri.^2);                    % the area of the cylinder 
I=((pi/4)*(r_og.^4-ri.^4));              % the second moment of inertia 
E=2*10^11;                               %young module
 %% formulation for the local K matrix then insert it in the global Matrix 
 kii=zeros(n1,n1);
 K=zeros(n1,n1);
 
 ki=(2*E*I/Li.^3)*[6 (3*Li) -6 (3*Li);(3*Li) (2*Li^2) (-3*Li) (Li.^2);-6 (-3*Li) 6 (-3*Li);(3*Li) (Li^2) (-3*Li) (2*Li.^2)];
 %%
 for i=1:e
   % creating the local matrix for each element
   kii((2*i)-1:(2*i)+2,(2*i)-1:(2*i)+2)=ki;
   % adding the element into the global matrix 
   K=K+kii;
   kii=zeros(n1,n1);
 end
 
%% formulation for the local M matrix then insert it in the global Matrix 
 mii=zeros(n1,n1);
 M=zeros(n1,n1);
 
mi=(rho*Li*A/420)*[156 (22*Li) 54 (-13*Li);(22*Li) (4*Li^2) (13*Li) (-3*Li.^2);54 (13*Li) 156 (-22*Li);(-13*Li) (-3*Li^2) (-22*Li) (4*Li.^2)];
 for i=1:e
    mii((2*i)-1:(2*i)+2,(2*i)-1:(2*i)+2)=mi;
    %%%% adding the local mass element into the global matrix 
    M=M+mii;
    mii=zeros(n1,n1);
  
 end
%%
%  Boundary conditions (cantilever beam) fixed at one end 
% the 1 st two node describing the movement of the lower end of the mast
% (translation and rotation ).
 
    K(1:2,:)=[]; 
    K(:,1:2)=[];
    M(1:2,:)=[];
    M(:,1:2)=[];
  
%% find the eigen value and the natural frequency for the cantilever beam 
 D=eig(K,M);
 omega=sqrt(D);
 f=omega/(2*pi);
 f_fund_num=f(1)
 %%
% the analytical solution 
f_fund_anly=((1.8751^2)/(2*pi*L^2))*sqrt((E*I)/(rho*A))

