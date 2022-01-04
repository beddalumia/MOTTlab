%% BSD 3-Clause License
% 
% Copyright (c) 2020, Gabriele Bellomia
% All rights reserved.

function I = LuttingerIntegral(w,Sigma_loc,gloc)
%% Computes the Luttinger sum-rule, as defined in PRB 102 081110
%  Input:
%       w           : real valued array, \omega domain
%       Sigma_loc   : complex valued array, \Sigma(\omega) function
%       gloc        : complex valued array, G_{loc}(\omega) function
%  Output:
%       I = \frac{2}{\pi}\Im\int_{-\infty}^{0}[ dw G_loc(w) d\Sigma(w)/dw ]
dw = w(2)-w(1);
dSigma = diff(Sigma_loc);
dSigma_Neg = dSigma(w<=0);
G_Neg = gloc(w<=0);
integrand = imag(G_Neg.*dSigma_Neg);
I = 2/pi*(dw*sum(integrand));
%plot(w(w<=0),integrand);
end

