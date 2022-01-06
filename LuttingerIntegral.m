%% BSD 3-Clause License
% 
% Copyright (c) 2020, Gabriele Bellomia
% All rights reserved.

function I = LuttingerIntegral(w,sloc,gloc)
%% Computes the Luttinger sum-rule, as defined in PRB 102 081110
%  Input:
%       w           : real valued array, \omega domain
%       sloc        : complex valued array, \Sigma(\omega) function
%       gloc        : complex valued array, G_{loc}(\omega) function
%  Output:
%       I = \frac{2}{\pi}\Im\int_{-\infty}^{0} dw G_loc(w) d\Sigma(w) / dw
%         = \frac{1}{\pi}\Im\int_{-\infty}^{\infty}dwG_loc(w)d\Sigma(w)/dw
dw  = w(2)-w(1);
ds  = diff(sloc);
g_  = gloc(1:end-1);
integrand  = smoothdata(imag(g_.*ds));
I = 1/pi*(dw*sum(integrand));
%plot(w_,integrand);
end

