function [A,B,C,D] = TM_cyl(r1,r2,L,f)
%% calculates the elements of a
%% Transmission matrix of a cylinder
%see e.g.
%Chaigne & Kergomard (2016): "Acoustics of Musical Instruments"
%Chap. 4.5; Eq.(4.28)
%https://link.springer.com/book/10.1007/978-1-4939-3679-3
%Timo Grothe, HfM Detmold, ETI 23.03.2023

global eta gamma Pr rho c

r = r1;

%characterisitic impedance
Zc = rho*c/(pi*r^2);

omega = 2*pi*f;% [rad/s]
k = omega/c;% [1/m]
% viscous boundary layer thickness
deltav = sqrt(eta./omega/rho); %[m]

% visco-thermal loss layer thickness
alphaprime = 1/2*sqrt(2)*((gamma - 1)/sqrt(Pr) + 1).*deltav;%[m]

req = r;%[m]

kc = k.*(1+alphaprime./req*(1-1i));

% %% Transmission matrix elements
A =        cos(kc*L);
B = 1i*Zc.*sin(kc*L);
C = 1i/Zc.*sin(kc*L);
D =        cos(kc*L);

