function [I,Lfig] = luttinger(w,sloc,gloc)
%% Computes the Luttinger sum-rule, as defined in PRB 90 075150
%
%  Input:
%       w           : real valued array, \omega domain
%       sloc        : complex valued array, \Sigma(\omega) function
%       gloc        : complex valued array, G_{loc}(\omega) function
%  Output:
%       I = \frac{1}{\pi}\Im\int_{-\infty}^{\infty}dwG_loc(w)d\Sigma(w)/dw
%
%% BSD 3-Clause License
% 
%  Copyright (c) 2022, Gabriele Bellomia
%  All rights reserved.
                                                               global DEBUG

% The original definition, as given by Logan, Tucker and Galpin, prescribes
% to integrate over the negative semiaxis:
%
%   I = \frac{2}{\pi} \Im\int_{-\infty}^{0} dw G_loc(w) d\Sigma(w) / dw
%      
% Nevertheless we find most useful to exploit particle-hole symmetry and
% thus integrate over the whole real-axis. Care has to be taken whenever
% the frequency domain contains the origin, for therein resides a nasty
% pole, arising from the derivative of the self-energy. Thus we divide 
% the integrand in (w < 0) and (w > 0) parts, and integrate the union.
% Failing to exclude the (w = 0) point, if exist, leads to a diverging
% uncontrolled result.

wl = w(w<0);        % Build w < 0 patch
gl = gloc(w<0);
sl = sloc(w<0); 
wl = wl(1:end-1);                   
gl = gl(1:end-1);
dl = diff(sl);

wr = w(w>0);        % Build w > 0 patch
gr = gloc(w>0); 
sr = sloc(w>0);
wr = wr(1:end-1); 
gr = gr(1:end-1);       
dr = diff(sr);    

w  = [wl,wr];       % Join the patches
ds = [dl,dr];
g  = [gl,gr];

integrand  = imag(g.*ds);
I = 1/pi*sum(integrand);
I = abs(I);
                                                                   if DEBUG
Lfig = figure("Name",'Luttinger integrand','Visible','off');
plot(w,integrand);
xlabel('\omega');
ylabel('Im[G(\omega)\partial\Sigma/\partial\omega]');
ylim([-0.1,0.2]);
                                                                   end
end


