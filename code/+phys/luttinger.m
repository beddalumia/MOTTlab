function [I,Lfig] = luttinger(w,sloc,gloc)
%% Computes the Luttinger sum-rule, as defined in PRB 102 081110
%
%  Input:
%       w           : real valued array, \omega domain
%       sloc        : complex valued array, \Sigma(\omega) function
%       gloc        : complex valued array, G_{loc}(\omega) function
%  Output:
%       I = \frac{2}{\pi}\Im\int_{-\infty}^{0} dw G_loc(w) d\Sigma(w) / dw
%         = \frac{1}{\pi}\Im\int_{-\infty}^{\infty}dwG_loc(w)d\Sigma(w)/dw
%
%% BSD 3-Clause License
% 
%  Copyright (c) 2020, Gabriele Bellomia
%  All rights reserved.
                                                               global DEBUG
                                                               w = w(w<=0);
                                                               sloc = sloc(w<=0);
                                                               gloc = gloc(w<=0);
ds  = diff(sloc); % dSigma -> already eliminates dw!
g_  = gloc(1:end-1);
integrand  = imag(g_.*ds);
I = 2/pi*sum(integrand);
I = abs(I); % J. Phys.: Condens. Matter 28 (2016) 025601
                                                                   if DEBUG
Lfig = figure("Name",'Luttinger integrand','Visible','off');
plot(w(1:end-1),integrand);
xlabel('\omega');
ylabel('Im[G(\omega)\partial\Sigma/\partial\omega]');
ylim([-0.1,0.2]);
                                                                   end
end


