%=============================================================================| 
% # Copyright (C) 2017 Dr.-Ing. Arun Raina (E-Mail: arunraina@icloud.com)
%
% This matlab script is part of the code used for the paper, 
% "Analysis of electro-permeation of hydrogen in metallic alloys".
% DOI: https://doi.org/10.1098/rsta.2016.0409
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This code is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with this code. If not, see <http://www.gnu.org/licenses/>.
%=============================================================================| 
function map_BDH_BN
clc; clear all; close all;
format long e

bmg = 40;
NTg1 = logspace(-5,-4,16);
NTg2 = logspace(-4,-3,13);
NTg3 = logspace(-3,-2,13);
NTg = unique([NTg1 NTg2 NTg3]);
% btaui(bDH_i,bN_i) : matrix of normalised time constant
for i = 1:bmg
    [btaui] = bDH_bN_datF(NTg(i),bmg);
    btau(:,i) = btaui; i
end
dlmwrite('filename.dat',btau,'delimiter',' ','precision','%5.4d')

% =====================================================
% MAIN FUNCTION
% =====================================================
function [btaui] = bDH_bN_datF(NTa,bmg)

% =====================================================
% define global quantities 
% =====================================================
global Q D0 R T0 l phi NL bet NT alp DH tL0 bm

% =====================================================
% List of input parameters for a ferritic steel
% =====================================================
Q   = 6680;             % Lattice enthalpy [J/mol]
D0  = 2.33e-7;          % Diffusion prefactor [m2/s]
R   = 8.3144;           % Gas constant [J/K/mol]
T0  = 293;              % Temperature [K]
m   = 0;                % PDE related plane slab
l   = 5e-3;             % Specimen thickness [m]
% tf  = 1e+6;           % Simulation final time [s]
phi = 0.0;              % heating rate [K/sec]
NL  = 8.46e+28;         % NILS density [atoms/m3]      
bet = 1;                % No. of NILS per lattice atom
NT  = NL*NTa;           % Traps site density [atoms/m3]
bm  = bmg; 
bDH = linspace(10,25,bm);
DH  =-R*T0*bDH;         % Traps binding energy [J/mol]    
alp = 1;                % No. of atoms per trap site   
tL0 = 1e-6;             % Initial lat. occpancy ratio

% =====================================================
% change total simulation time bQ=2.75
% =====================================================
if NTa >=9e-6 && NTa <2e-5
    tf=3e4;
elseif NTa >=2e-5 && NTa <8e-5
    tf=1e5;
elseif NTa >=8e-5 && NTa <8e-4
    tf=1e6;    
elseif NTa >=8e-4 && NTa <=1e-2
    tf=1e7;     
end

% =====================================================
% Space & time discretization
% =====================================================
nx = 1e2;
nt = 1e4;
xa = linspace(-l/2,l/2 ,nx);
ta = linspace(0,tf,nt);

% =====================================================
% Normalised space and time
% =====================================================
x  = xa./l;
t  = (D0/l^2).*ta;
dx = x(2)-x(1);

% =====================================================
% Solving ODE for \bar\theta_L 
% =====================================================
sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);
u = zeros(nt,nx,bm);
for i = 1:bm
    u(:,:,i)=sol(:,:,i);  
end

% =====================================================
% Computing normalised Flux at x = +l/2 at each time step
% =====================================================
J=zeros(nt,bm);
for j = 1:bm
    for i = 1:nt
        J(i,j) = -(u(i,nx,j)-u(i,nx-1,j))/dx;
    end
end

% =====================================================
% Find normalised time constant when J=0.632
% =====================================================
btaui(1:bm)=0;
diff(1)=0.632;
for j = 1:bm
    for i = 2:nt
        diff(i) = abs(J(i,j)-0.632);
        if diff(i)>diff(i-1)
            btaui(j)=t(i-1);
            jxxx(j)=J(i-1,j);
            break
        end
    end
end

% =====================================================
% Plotting
% =====================================================
% figure('units','normalized','position',[.1 .1 .4 .25])
% for j = 1:bm
%     semilogx(t,J(:,j)); hold on
%     semilogx(btaui(j),jxxx(j),'xk')
% end
% ylabel('J')
% xlabel('T')
% grid on
% set(gca,'Fontsize',13)
clear global Q D0 R T0 l phi NL bet NT alp DH tL0

% =====================================================
% PDE coefficients
% =====================================================
function [c,f,s] = pdefun(x,t,u,DuDx)
global Q D0 R T0 l phi NL bet NT alp DH tL0 bm

% =====================================================
% Normalised quantities
% =====================================================
bQ   = Q/(R*T0);
bphi = phi*l^2/(T0*D0);
bT   = 1 + bphi*t;
bN   = (alp*NT)/(bet*NL);
bDH  = DH./(R*T0);
K    = exp(-bDH./bT); 
bDL  = exp(-bQ/bT);
% 
for i = 1:bm
    Di(i) = 1 + (bN*K(i))/(1+tL0*K(i)*u(i))^2;
    St(i) = -(bphi*u(i))/(bT^2)*(bN*K(i)*bDH(i))/(1+tL0*K(i)*u(i))^2;
end
%
c = Di(1:bm)'; 
f = repmat(bDL,1,bm)'.*DuDx;
s = St(1:bm)';

% =====================================================
% Initial condition
% =====================================================
function u0 = icfun(x)
global bm
u0 = zeros(1,bm)';

% =====================================================
% Boundary conditions
% =====================================================
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global bm
pl = ul(1:bm)-1.0;
ql = zeros(1,bm)';
pr = ur(1:bm); 
qr = zeros(1,bm)';
