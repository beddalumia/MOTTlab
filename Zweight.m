%% BSD 3-Clause License
% 
% Copyright (c) 2020, Gabriele Bellomia
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

