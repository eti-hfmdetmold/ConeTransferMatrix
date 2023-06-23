function [A,B,C,D] = TM_con_nederveen(r1,r2,L,f)
%% calculates the elements of a
%% Transmission matrix of a cone
%Ref: Grothe, Baumgart, Nederveen (2023): "A Transfer Matrix for the Input Impedance of 
%weakly tapered Cones as of Wind Instruments" 
%https://arxiv.org/abs/2303.127505 (preprint)
%Timo Grothe, HfM Detmold, ETI 23.03.2023

global eta gamma Pr rho c

if r1 == r2
    r2 = r1+eps;
end

%characteristic impedance
Zc1 = rho*c/(pi*r1^2);

%taper:
m = (r2-r1)/L;%[-]
%input distance from apex
x1 = r1/m;%[m]
%output distance from apex
x2 = r2/m;%[m]

omega = 2*pi*f;% [rad/s]
k = omega/c;% [1/m]
% viscous boundary layer thickness 
deltav = sqrt(eta./omega/rho); %[m]

% visco-thermal loss layer thickness
alphaprime = 1/2*sqrt(2)*((gamma - 1)/sqrt(Pr) + 1).*deltav;%[m]

%%Nederveens helper values f and g (with Ci and si as defined in Abramovitz
%%and Stegun.

% option (1)
% uses Matlabs built-in functions sinint.m and cosint.m
% fAS =@(x)  cosint(2*k*x).*sin(2*k*x) - (sinint(2*k*x)-pi/2).*cos(2*k*x); %[-]
% gAS =@(x) -cosint(2*k*x).*cos(2*k*x) - (sinint(2*k*x)-pi/2).*sin(2*k*x); %[-]
% 
% f1 = fAS(x1);f2 = fAS(x2);
% g1 = gAS(x1);g2 = gAS(x2);

% option (2) 
% Uses approximations of GalSim: The modular galaxy image simulation
% toolkit (Rowe et al. Astronomy and Computing, 2015), which prove to be
% much faster than the Matlab built-ins cosint and sinint

[f1,g1] = f_g_approximations(2*k*x1);
[f2,g2] = f_g_approximations(2*k*x2);

% loss corrections (complex)
ag1 = 1/2*sqrt(2)*deltav.*(1 + 2.*g1.*(k*x1).^2*((gamma - 1)/sqrt(Pr) + 1))./r1*(1i-1);
ag2 = 1/2*sqrt(2)*deltav.*(1 + 2.*g2.*(k*x2).^2*((gamma - 1)/sqrt(Pr) + 1))./r2*(1i-1);

af1 = 1/2*sqrt(2)*deltav.*(2 - 2*f1.*(k*x1).*((gamma - 1)/sqrt(Pr) + 1))./r1*(1i-1);
af2 = 1/2*sqrt(2)*deltav.*(2 - 2*f2.*(k*x2).*((gamma - 1)/sqrt(Pr) + 1))./r2*(1i-1);


%effective cone radius (frequency dependent)
%reff = (r2-r1)./(g2-g1+log(r2/r1));%[m] 
reff = (r2-r1)./(g2-g1+log1p((r2-r1)./r1));%[m] %use log1p (thanks Augustin)

%argument of sin , cos
sigma = (k.*L.*(1+alphaprime./reff*(1-1i)));

% %% Transmission matrix elements
A =     r2/r1*cos(sigma).*(1 + af2)...
    - 1./(k.*x1).*sin(sigma).*(1 + ag2);

B =  Zc1.*(r1/r2.*sin(sigma)...
    ).*1i;

C = 1./Zc1.*(...
    r2/r1*sin(sigma).*(1 + af1).*(1 + af2)...
    + 1./(k.*x1).*(...
    1./(k.*x1).*sin(sigma).*(1 + ag1).*(1 + ag2)...
    - r2/r1.*cos(sigma).*(1 + af2).*(1 + ag1)...
    + cos(sigma).*(1 + af1).*(1 + ag2)...
    )...
    ).*1i; 


D = (r1/r2).*(              cos(sigma).*(1 + af1) +...
    1./(k*x1).*sin(sigma).*(1 + ag1) ...
    );
end

function [f,g] = f_g_approximations(x)
%initialize
f = zeros(size(x));
g = zeros(size(x));

%Rowe2015 https://arxiv.org/abs/1407.7676

idxsmall = find(x<4);
idxlarge = find(x>=4);

xsmall = x(idxsmall);
xlarge = x(idxlarge);

x = xsmall;
%Eq. (B.5)
Si = x.*...
    ( 1 -4.54393409816329991e-2 * x.^2 + 1.15457225751016682e-3 * x.^4 - 1.41018536821330254e-5 * x.^6 + ...
         9.43280809438713025e-8 * x.^8 - 3.53201978997168357e-10 * x.^10 + 7.08240282274875911e-13 * x.^12 -...
         6.05338212010422477e-16 * x.^14 )./ ...
    ( 1 + 1.01162145739225565e-2 * x.^2 + 4.99175116169755106e-5 * x.^4 + 1.55654986308745614e-7 * x.^6 + ...
         3.28067571055789734e-10 * x.^8 + 4.5049097575386581e-13 * x.^10 + 3.21107051193712168e-16 * x.^12);

gamma = 0.57721566490153286060651209008240243104215933593992;
% Eq. (B.6)
Ci = gamma + log(x) + x.^2 .* ...
    (-0.25 + 7.51851524438898291e-3 * x.^2 - 1.27528342240267686e-4  * x.^4  + 1.05297363846239184e-6 * x.^6 -...
             4.68889508144848019e-9 * x.^8 + 1.06480802891189243e-11 * x.^10 - 9.93728488857585407e-15 * x.^12)./...
    ( 1   + 1.1592605689110735e-2 * x.^2 + 6.72126800814254432e-5 * x.^4 + 2.55533277086129636e-7 * x.^6 +...
            6.97071295760958946e-10* x.^8 + 1.38536352772778619e-12* x.^10 + 1.89106054713059759e-15 * x.^12 +...
            1.39759616731376855e-18* x.^14);
        
% Eq. (B.7)
f(idxsmall) =  Ci.*sin(x)+(pi/2-Si).*cos(x);
% Eq. (B.8)
g(idxsmall) = -Ci.*cos(x)+(pi/2-Si).*sin(x);


x = xlarge;    
% Eq. (B.11)   
f(idxlarge) = 1./x .* ...
    (1 + 7.44437068161936700618e2 * x.^-2 + 1.96396372895146869801e5 * x.^-4 + 2.37750310125431834034e7 * x.^-6 +...
         1.43073403821274636888e9 * x.^-8 + 4.33736238870432522765e10 * x.^-10 + 6.40533830574022022911e11 * x.^-12 +...
         4.20968180571076940208e12 * x.^-14 + 1.00795182980368574617e13 * x.^-16 + 4.94816688199951963482e12 * x.^-18 -...
         4.94701168645415959931e11 * x.^-20)./...
    (1 + 7.46437068161927678031e2 * x.^-2 + 1.97865247031583951450e5 * x.^-4 + 2.41535670165126845144e7 * x.^-6 +...
         1.47478952192985464958e9 * x.^-8 + 4.58595115847765779830e10 * x.^-10 + 7.08501308149515401563e11 * x.^-12 +...
         5.06084464593475076774e12 * x.^-14 + 1.43468549171581016479e13 * x.^-16 + 1.11535493509914254097e13 * x.^-18);

     %Eq. (B.12)
g(idxlarge) = x.^-2 .*...
    (1 + 8.1359520115168615e2 * x.^-2 + 2.35239181626478200e5 * x.^-4 +3.12557570795778731e7 * x.^-6 +...
         2.06297595146763354e9 * x.^-8 + 6.83052205423625007e10 * x.^-10 + 1.09049528450362786e12 * x.^-12 +...
         7.57664583257834349e12 * x.^-14 + 1.81004487464664575e13 * x.^-16 + 6.43291613143049485e12 * x.^-18 - ...
         1.36517137670871689e12 * x.^-20)./ ...
   (1 + 8.19595201151451564e2 * x.^-2 + 2.40036752835578777e5 * x.^-4  + 3.26026661647090822e7 * x.^-6 +...
        2.23355543278099360e9 * x.^-8 + 7.87465017341829930e10 * x.^-10 + 1.39866710696414565e12 * x.^-12 + ...
        1.17164723371736605e13 * x.^-14 + 4.01839087307656620e13 * x.^-16 + 3.99653257887490811e13 * x.^-18);
    
end

