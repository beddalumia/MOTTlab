function im_sloc = sopt(A,F,U)
%% SOPT stands for Second Order Perturbation Theory:
%  it computes the imag of the 2nd order diagram for the self-energy
%
%% Input  
%
%  A  : Noninteracting spectral function, A(w) -> real array
%  F  : Fermi-Dirac distribution, f(w) -> real array
%  U  : Hubbard interaction, U -> real scalar
%
%% Output
%
%  iS : Im(Sigma(w)) -> real array
%
%% Theoretical Background at:
%
%  http://www.physics.rutgers.edu/~haule/681/Perturbation.pdf
%
%% BSD 3-Clause License
%
%  Copyright (c) 2020, Gabriele Bellomia
%  All rights reserved.

global DEBUG FAST
                                                          
  %% Particle distribution: A^+ = A*f [Hole distribution: A^- = A(1-f)]
     A = A.*F; 
     
if FAST 
    
  %% Optimized FFTW-based convolution
     P = math.fconv(A,A,'same');  % Polarization BUBBLE
     D = math.fconv(A,P,'same');  % 2nd-order DIAGRAM
     
else
    
  %% Built-in vectorized convolution
     P = conv(A,A,'same');        % Polarization BUBBLE       
     D = conv(A,P,'same');        % 2nd-order DIAGRAM   
     
end

  %% Enforce ph & half-filling to retrieve A^- contribution (holes)
     D = D + flip(D); % flip{v(1:end)}=v(end:1)
     
  %% Imaginary part of the Self-Energy according to SOPT
     im_sloc = -pi * U^2 * D;
    
if FAST && DEBUG
    
  %% Cross-Check
     P_test = conv(A,A,'same');
     err1 = abs(norm(P_test-P)/norm(P_test));
     if err1 > 10*eps %-> one order above machine precision
        fprintf(2,'Error on polarization bubble: %.16f \n',err1);
     end
     D_test = conv(A,P_test,'same');
     err2 = abs(norm(D_test-D)/norm(D_test));
     if err2 > 10*eps %-> idem
        fprintf(2,'Error on SOPT diagram: %.16f \n',err2);
     end
     
end


