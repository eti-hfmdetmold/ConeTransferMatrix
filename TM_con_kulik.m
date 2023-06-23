function [A,B,C,D] = TM_con_kulik(r1,r2,L,f)
%% calculates the elements of a
%% Transmission matrix of a cone
%Ref: Kulik (JASA-EL 2007): "Transfer matrix of conical waveguides with any geometric parameters for increased
%precision in computer modeling" 
%https://asa.scitation.org/doi/full/10.1121/1.2794865)
%Timo Grothe, HfM Detmold, ETI 23.03.2023

global eta gamma Pr rho c

%taper:
m = (r2-r1)/L;%[-]
%input distance from apex
x1 = r1/m;%[m]
%output distance from apex
x2 = r2/m;%[m]

L12 = L;

omega = 2*pi*f;% [rad/s]
k = omega/c;% [1/m]
% viscous boundary layer thickness
deltav = sqrt(eta./omega/rho); %[m]

% visco-thermal loss layer thickness
alphaprime = 1/2*sqrt(2)*((gamma - 1)/sqrt(Pr) + 1).*deltav;%[m]

%reff = (r2-r1)./(log(r2/r1));%[m] 
req = (r2-r1)./(log1p((r2-r1)./r1));%[m]

k1 = k.*(1+alphaprime./r1*(1-1i));
k2 = k.*(1+alphaprime./r2*(1-1i));
k12 = k.*(1+alphaprime./req*(1-1i));

theta1 = atan(k1*x1); theta2 = atan(k2*x2);
t1     = 1./sin(theta1); t2 = 1./sin(theta2);

% %% Transmission matrix elements
A =         t2.*cos(k12*L12-(theta2-pi/2));
B = 1i        .*sin(k12*L12);
C = 1i.*t1.*t2.*sin(k12*L12-theta2+theta1);
D =         t1.*cos(k12*L12+(theta1-pi/2));

