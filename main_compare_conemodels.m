clear all
close all
global eta gamma Pr rho c

%This script computes the input impedance of a weakly tapered, dissipative conical frustum
%using the transmission line method with three different transmission matrices.
%The results are compared to a numerical reference solution computed with
%1D finite-elements, computed with the python library openWind.
%
%Required functions:
% cons_ow.m
% TM_con_kulik.m
% TM_con_nederveen.m
% TM_cyl.m
% TML.m
%Reference data:
% closed_cone(ow).txt (from closed_cone.py)
%
% Timo Grothe, HfM Detmold, 22.06.2023

%% temperature
tc = 20;%[°C]

%% reference curve
%0) numerical solution from openWind https://openwind.inria.fr/ 
%(courtesy of Augustin Ernoult, computed @tc = 20°C)

data = load('closed_cone(ow).txt');% this result file is computed with open_cone.py 
f = data(:,1);
Zref = data(:,2)+1i*data(:,3);

%% physical constants (Ref. Chaigne&Kergomard Eq. 5.142 (from openwind))
[c, rho, delta, epsil, eta, gamma, kappa, Cp, nu] = cons_ow(tc);
Pr = nu^2;%[-]

%% geometry
%input, output radii, length L [m]
r1 = 2.1e-3;% 
r2 = 23.5e-3;
L = 3;


%% impedance at the ducts far end 
Z2 = 1e32; %[Pa*s/m^3] (closed end)

%% characteristic impedances
Zc1 = rho*c/(pi*r1^2);%[Pa*s/m^3]
Zc2 = rho*c/(pi*r2^2);%[Pa*s/m^3]

%% Transmission matrices
%1) Grothe, Baumgart, Nederveen (manuscript in preparation)
[A,B,C,D] = TM_con_nederveen(r1,r2,L,f);
Z_nederveen = (A.*Z2+B)./(C.*Z2+D);
%2) Kulik (Transfer matrix of conical waveguides with any geometric parameters for increased
%precision in computer modeling, JASA 2007)
[A,B,C,D] = TM_con_kulik(r1,r2,L,f);
Z_kulik = Zc1*(A.*Z2/Zc2+B)./(C.*Z2/Zc2+D); %Kuliks original transfer matrix is normalized
%3) analytical (e.g. Chaigne/ Kergomard: Acoustics of Musical Instruments, Springer 2016)
n = 6000;% number of cylindrical slices approximatig the cone
[A,B,C,D] = TML(r1,r2,L,f,n,'cyl','classical');
Z_cyls = (A.*Z2+B)./(C.*Z2+D);


%% plot
subplot(2,2,1)
semilogy(f,real(Zref))
hold on
plot(f,real(Z_nederveen))
plot(f,real(Z_kulik))
plot(f,real(Z_cyls))

subplot(2,2,2)
plot(f,imag(Zref))
hold on
plot(f,imag(Z_nederveen))
plot(f,imag(Z_kulik))
plot(f,imag(Z_cyls))

subplot(2,2,3)
ax = gca; ax.ColorOrderIndex = 2;
hold on
semilogy(f,abs((real(Z_nederveen)-real(Zref))./(real(Zref))),'.')
semilogy(f,abs((real(Z_kulik)-real(Zref))./(real(Zref))),'.')
semilogy(f,abs((real(Z_cyls)-real(Zref))./(real(Zref))),'.')
set(gca,'yscale','log')

subplot(2,2,4)
ax = gca; ax.ColorOrderIndex = 2;
hold on
semilogy(f,abs((imag(Z_nederveen)-imag(Zref))./(imag(Zref))),'.')
semilogy(f,abs((imag(Z_kulik)-imag(Zref))./(imag(Zref))),'.')
semilogy(f,abs((imag(Z_cyls)-imag(Zref))./(imag(Zref))),'.')
set(gca,'yscale','log')

subplot(2,2,1)
legend('1D FEM (ow)','TM single cone (nederveen)','TM single cone (kulik)','TM cyl.slices',...
    'location','southeast')

subplot(2,2,1)
ylabel('real (Z) [Pa*s/m^3]')
xlabel('f [Hz]')
subplot(2,2,2)
ylabel('imag (Z) [Pa*s/m^3]')
xlabel('f [Hz]')
subplot(2,2,3)
ylabel('real(Z-Z_{ow})/real(Z_{ow}) [-]')
xlabel('f [Hz]')
subplot(2,2,4)
ylabel('imag(Z-Z_{ow})/imag(Z_{ow}) [-]')
xlabel('f [Hz]')