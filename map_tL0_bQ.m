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
function map_tL0_bQ
clc; clear all; close all;
format long e

bm = 40;
Q  = linspace(6.68e3,33.5e3,bm);
tL01 = logspace(-7,-6,16);
tL02 = logspace(-6,-5,13);
tL03 = logspace(-5,-4,13);
tL0  = unique([tL01 tL02 tL03]);
% Jssi(tL0_i,Q_j) : matrix of steady state flux
for i = 1:bm
    for j = 1:bm
        [Jss1, btau1] = bDH_bN_datF(tL0(i),Q(j));
        Jss(i,j)  = Jss1; i
        btau(i,j) = btau1; j
    end
end
dlmwrite('valuesJ.dat',Jss,'delimiter',' ','precision','%5.4d')
dlmwrite('valuest.dat',btau,'delimiter',' ','precision','%5.4d')

% =====================================================
% MAIN FUNCTION
% =====================================================
function [Jssi, btaui] = bDH_bN_datF(tL0i,Qj)

% =====================================================
% define global quantities 
% =====================================================
global Q D0 R T0 l phi NL bet NT alp DH tL0

% =====================================================
% List of input parameters
% =====================================================
Q   = Qj;                % Lattice enthalpy [J/mol]
D0  = 2.33e-7;           % Diffusion prefactor [m2/s]
R   = 8.3144;            % Gas constant [J/K/mol]
T0  = 293;               % Temperature [K]
m   = 0;                 % PDE related plane slab
l   = 5e-3;              % Specimen thickness [m]
% tf  = 1e+4;            % Simulation final time [s]
phi = 0.0;               % heating rate [K/sec]
NL  = 8.46e+28;          % NILS density [atoms/m3]      
bet = 1;                 % No. of NILS per lattice atom
NT  = NL*1e-4;           % Trap site density [atoms/m3]
bDH = 15;                % <<--   
DH  =-R*T0*bDH;          % Trap binding energy [J/mol]    
alp = 1;                 % No. of atoms per trap site   
tL0 = tL0i;

% =====================================================
% change total simulation time (bDH=15)
% =====================================================
if Qj >=6.6e3 && Qj <10e3
    tf=1e6;
elseif Qj >=10e3 && Qj <14e3
    tf=1e7;
elseif Qj >=14e3 && Qj <20e3
    tf=1e8;    
elseif Qj >=20e3 && Qj <24e3
    tf=1e9;   
elseif Qj >=24e3 && Qj <30e3
    tf=1e10;       
elseif Qj >=30e3
    tf=1e11;         
end

% =====================================================
% Space & time discretization
% =====================================================
nx = 1e2;
nt = 1e6; % for bDH=10 use 1e4 
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
u(:,:)=sol(:,:);  

% =====================================================
% Computing normalised Flux at x = +l/2 at each time step
% =====================================================
for i = 1:nt
    J(i) = -(u(i,nx)-u(i,nx-1))/dx;
end

% =====================================================
% Find normalised time constant when J=0.632
% =====================================================
diff(1)=0.632;
for i = 2:nt
    diff(i) = abs(J(i)-0.632);
    if diff(i)>diff(i-1)
        btaui=t(i-1);
        jxxx =J(i-1);
        break
    end
end
% steady state flux at x = +l/2 in [mol/m^2/s]
Jssi=D0*exp(-Q/R/T0)*bet*NL*tL0/l/6.02214e+23;

% =====================================================
% Plotting
% =====================================================
% figure('units','normalized','position',[.1 .1 .4 .25])
% semilogx(t,J); hold on
% semilogx(btaui,jxxx,'xk')
% ylabel('J')
% xlabel('T')
% grid on
% set(gca,'Fontsize',13)
clear global Q D0 R T0 l phi NL bet NT alp DH tL0

% =====================================================
% PDE coefficients
% =====================================================
function [c,f,s] = pdefun(x,t,u,DuDx)
global Q D0 R T0 l phi NL bet NT alp DH tL0

% =====================================================
% Normalised quantities
% =====================================================
bQ   = Q/(R*T0);
bphi = phi*l^2/(T0*D0);
bT   = 1 + bphi*t;
bN   = (alp*NT)/(bet*NL);
bDH  = DH/(R*T0);
K    = exp(-bDH/bT);
bDL  = exp(-bQ/bT);
% 
Di = 1 + (bN*K)/(1+tL0*K*u)^2;
St = -(bphi*u)/(bT^2)*(bN*K*bDH)/(1+tL0*K*u)^2;
%
c = Di; 
f = bDL*DuDx;
s = St;

% =====================================================
% Initial condition
% =====================================================
function u0 = icfun(x)
u0 = 0;

% =====================================================
% Boundary conditions
% =====================================================
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
pl = ul-1.0;
ql = 0;
pr = ur; 
qr = 0;
