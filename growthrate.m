function L = growthrate(a1, G1, a2, G2, b, Re, k, model)
%growthrate Summary of this function goes here
%   Detailed explanation goes here
if G1 < abs(G2)
    print('\Gamma_1 should greater than |\Gamma_2|\n');
    return
end
eta1 = Fetai(a1, b, Re, k);
eta2 = Fetai(a2, b, Re, k);
omega1 = Fomegai(a1, b, k, model);
omega2 = Fomegai(a2, b, k, model);
psi = Fpsi(b, k);
Chi = FChi(b, k);
Lambda = G2 / G1;
L = [eta1              , 1+omega1, 0                   , Lambda*psi;
     -1-2*Lambda-omega1, eta1    , Lambda*Chi          , 0;
     0                 , psi     , eta2                , Lambda*(1+omega2);
     Chi               , 0       , -2-Lambda*(1+omega2), eta2];
end

function res = Fetai(a, b, Re, k)
ka = k * a;
res = - (2*pi/Re) * (a/b)^(-2) * (1.54*ka+ka*ka);
end

function res = Fomegai(a, b, k, model)
gamma = 0.577;
if model==0 % Rankin vortex, uniform
    K = 0.25
elseif model==1 %Lamb-Oseen vortex
    K = -0.058;
end
ka = k * a;
kb = k * b;
%leweke 2016
term1 = kb * kb / (2 + 0.995 * ka + 0.438 * ka * ka);
term2 = log((2 + 2.151 * ka) / ka) + K - gamma;
res = term1 * term2;
% Widnall 1971
% res = kb * kb / 2 * (log(2/ka) + K - gamma);
end

function res = Fpsi(b, k)
kb = b * k;
res = kb * kb * besselk(0, kb) + kb * besselk(1, kb);
end

function res = FChi(b, k)
kb = k * b;
res = kb * besselk(1, kb);
end

