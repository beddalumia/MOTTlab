%% BSD 3-Clause License
% 
% Copyright (c) 2020, Gabriele Bellomia
% All rights reserved.

function Z = Zweight(w,Sigma_loc)
%% ZWEIGHT extracts the Landau quasiparticle weight from \Sigma(\omega)
%  Input:
%       w           : real valued array, \omega domain
%       Sigma_loc   : complex valued array, \Sigma(\omega) function
w_th = 0.01; % Threshold to define proper fitting domain! (to be improved)
xToFit = w(abs(w)<=w_th);
yToFit = real(Sigma_loc(abs(w)<=w_th));
LinearModel = polyfit(xToFit,yToFit,1);
%% Uncomment to debug / fine-tuning of w_th %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yModel = LinearModel(1)*w + LinearModel(2);                             %
% figure("Name",'Debug on Z determination');                              %
% plot(w,real(Sigma_loc),':'); hold on                                    %
% plot(w, yModel,'r');                                                    %
% ylim([min(real(Sigma_loc)),max(real(Sigma_loc))]);                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nozieres Theorem: Re{\Sigma(w)} = A + (1-1/Z)w + O(w^3)
Z = 1/(1-LinearModel(1)); 
end

