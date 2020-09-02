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

function Im_2ndDiagram = PT(A0,f,U)
%% PT stands for Perturbation Theory:
%  it computes the imaginary part of the 2nd order diagram (bubble)
%% Input  
%  A0: Noninteracting spectral function, A_0(\omega) -> real-value array
%  f: Fermi-Dirac distribution, f(\omega) -> real-value array
%  U: On-site interaction, U -> real float value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NB. We take advantage of half-filling condition and discard A^-(\omega)
Ap = A0.*f; % Occupied States Distribution: A^+(\omega)
App = conv(Ap, Ap, 'same');     % Fast built-in convolution
Appp = conv(Ap, App, 'same');   % Fast built-in convolution
pppA = flip(Appp); % flip([1 2 3]) == [3 2 1]
Im_2ndDiagram = -pi*U^2*(Appp + pppA);
%% Theoretical Background at:
%  http://www.physics.rutgers.edu/~haule/681/Perturbation.pdf
end

