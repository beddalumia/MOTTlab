function [Z,Zfig] = zetaweight(w,sloc)
%% ZWEIGHT extracts the Landau quasiparticle weight from \Sigma(\omega)
%
%  Input:
%       w    : real valued array, \omega domain
%       sloc : complex valued array, \Sigma(\omega) function
%  Output:
%       Z    : 1 / ( 1 - d/dw[Re\Sigma]@w=0 )
%
%% BSD 3-Clause License
% 
%  Copyright (c) 2020, Gabriele Bellomia
%  All rights reserved.

                                                               global DEBUG

w_th = 0.1; % Threshold to define proper fitting domain! (to be improved)
xToFit = w(abs(w)<=w_th);
yToFit = real(sloc(abs(w)<=w_th));
LinearModel = polyfit(xToFit,yToFit,1);

                                                                   if DEBUG
%% debug / fine-tuning of w_th
yModel = LinearModel(1)*w + LinearModel(2);                             
Zfig = figure("Name",'Debug on Z determination','Visible','off');               
plot(w,real(sloc),':'); hold on                                         
plot(w, yModel,'r');
xlabel('\omega');
legend('Re\Sigma(\omega)','d/dw[Re\Sigma]@\omega=0')
a = unique(min(real(sloc)));
b = unique(max(real(sloc)));
if a<b, ylim([a,b]); end                                
                                                                   end

%% Nozieres Theorem: Re{\Sigma(w)} = A + (1-1/Z)w + O(w^3)
Z = 1/(1-LinearModel(1)); 
end


