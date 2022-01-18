function Im_2ndDiagram = SOPT(A0,f,U)
%% SOPT stands for Second Order Perturbation Theory:
%  it computes the imaginary part of the 2nd order diagram (bubble)
%
%% Input  
%
%  A0 : Noninteracting spectral function, A_0(\omega) -> real-value array
%  f  : Fermi-Dirac distribution, f(\omega) -> real-value array
%  U  : On-site interaction, U -> real scalar
%
%% Theoretical Background at:
%
%  http://www.physics.rutgers.edu/~haule/681/Perturbation.pdf
%
%% BSD 3-Clause License
%
%  Copyright (c) 2020, Gabriele Bellomia
%  All rights reserved.

%% NB. We take advantage of half-filling condition and discard A^-(\omega)
Ap = A0.*f; % Occupied States Distribution: A^+(\omega)
App = conv(Ap, Ap, 'same');     % Fast built-in convolution
Appp = conv(Ap, App, 'same');   % Fast built-in convolution
pppA = flip(Appp); % flip([1 2 3]) == [3 2 1]
Im_2ndDiagram = -pi*U^2*(Appp + pppA);
end

