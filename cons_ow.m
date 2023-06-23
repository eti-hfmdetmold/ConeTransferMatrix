function [c, rho, delta, epsil, eta, gamma, kappa, Cp, nu] = cons_ow(tc)
%% INPUT
%tc :temperature [�C]
%% OUTPUT
%air and wave propagation constants [c rho delta epsil]
%          c      :speed of sound in free space [m/s]
%          rho    :air density [kg/m^3]
%          delta  :constant for thermal and viscous duct losses [m/sqrt(s)]
%          epsil  :constant for thermal and viscous duct losses [sqrt(s)]

% Augustin Ernoult 09.2022


% Ref Chaigne & Kergo, (Chap. 5, p.241).

cal2joule = 4.184 ;% conversion cal to joule

T0 = 273.15; % �C - Kelvin conversion

% Constant imposed
Cp_cal = 240; % Cal/(kg.°C)  specific heat with constant pressure
gamma = 1.402; %heat capacity ratio (adiabatic index) Cp/Cv

T     = tc + T0; % temp in K

rho   =  1.2929 * T0 / T; % air density in kg/m3
c     =  331.45 * sqrt(T / T0); % air sound speed in m/s
eta    = 1.708e-5 * (1 + 0.0029 * tc); % viscosity in kg/(m.s) or [N*s/m²]
kappa_cal = 5.77 * 1e-3 * (1 + 0.0033 * tc); %thermal conductivity in Cal/(m s °C)
khi   = 1 / (rho * c ^ 2);  % isentropic compressibility in kg/(m.s)
c_lt  = kappa_cal / (rho * Cp_cal); % step of computation
c_lv  = eta / rho; % step of computation
lt    = c_lt / c; %characteristic width of thermal effect in m
lv    = c_lv / c; % characteristic width of viscous effect in m

Pr = lv/lt; % Prandtl number = eta*Cp/kappa
nu = sqrt(Pr); % square root of Prandtl number
Cp = Cp_cal*cal2joule; % J/(K.°C)
kappa = kappa_cal * cal2joule; % J/(m.s.°C)

%%
% reference
% Benade1968 JASA
%with Eq(4) and (5)
delta = 1./2*sqrt(1./pi*eta./rho).*(1 +(gamma-1)./nu); %Eq 17 (c) 
epsil = 2*pi./c.*delta; %Eq 17 (d)



