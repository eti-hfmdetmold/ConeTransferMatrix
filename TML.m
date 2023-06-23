function [A,B,C,D] = TML(r1,r2,L,f,n,type,method)
%Transmission-Line Scheme
%Timo Grothe, HfM Detmold, ETI
%23.03.2023

xi = linspace(0,L,(n+1));
ri = interp1([0,L],[r1,r2],xi);
%initialize for frequencies
nfreq = size(f,1);
%initialize transfer matrix
T = zeros(2,2,nfreq);
T(1,1,:) = 1;
T(2,2,:) = 1;

%% build transmission line
for j = 1:n
    if strcmp(type,'cyl')
        %% cylindrical slices
        rj1 = (ri(j)+ri(j+1))/2;
        rj2 = rj1;
    elseif strcmp(type,'con')
        %% conical slices
        rj1 = ri(j);
        rj2 = ri(j+1);
    end
    %length of slice
    Lj  = xi(j+1)-xi(j);
    %% Transmission Matrix elements
    if     strcmp(method,'classical')
        [Aj,Bj,Cj,Dj] = TM_cyl(rj1,rj2,Lj,f);
    elseif strcmp(method,'nederveen')
        [Aj,Bj,Cj,Dj] = TM_con_nederveen(rj1,rj2,Lj,f);
    elseif strcmp(method,'nederveen_') %normalized version
        [Aj,Bj,Cj,Dj] = TM_con_nederveen_(rj1,rj2,Lj,f);
    elseif strcmp(method,'kulik')
        [Aj,Bj,Cj,Dj] = TM_con_kulik(rj1,rj2,Lj,f);
    end
    %transmission matrix of current segment
    Tj(1,1,:) = reshape(Aj,1,1,[]);
    Tj(1,2,:) = reshape(Bj,1,1,[]);
    Tj(2,1,:) = reshape(Cj,1,1,[]);
    Tj(2,2,:) = reshape(Dj,1,1,[]);
    %update transmission matrix of chain by transmission matrix of current segment
    for i = 1:nfreq
        T(:,:,i) = T(:,:,i)*Tj(:,:,i);
    end
end
%% resulting TM coefficients
A = squeeze(T(1,1,:));
B = squeeze(T(1,2,:));
C = squeeze(T(2,1,:));
D = squeeze(T(2,2,:));