function [Im_2ndDiagram, dbf1, dbf2] = SOPT(A0,f,U)
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
                                                          global DEBUG FAST
                                                          
    %% We take advantage of half-filling condition and discard A^-(\omega)
    Ap   = A0.*f;              % Occupied States Distribution: A^+(\omega)
  if FAST
    App  = math.fconv(Ap,Ap,'same');  % Optimized FFTW-based convolution
    Appp = math.fconv(Ap,App,'same'); % Optimized FFTW-based convolution
  else
    App  = conv(Ap,Ap,'same');        % Built-in vectorized convolution
    Appp = conv(Ap,App,'same');       % Built-in vectorized convolution      
  end
    pppA = flip(Appp);  % flip([1 2 3]) == [3 2 1]
    Im_2ndDiagram = -pi*U^2*(Appp + pppA);
  if FAST && DEBUG
    App_test = conv(Ap,Ap,'same');
    dbf1 = figure("Name",'dbf1','Visible','off');
    plot(length(Ap)/2:length(Ap),App(length(Ap)/2:length(Ap))); hold on; 
    plot(length(Ap)/2:length(Ap),App_test(length(Ap)/2:length(Ap)),'--');
    err1 = abs(norm(App_test-App)/norm(App_test));
    if err1 > 10*eps %-> one order above machine precision
       fprintf('Error on 1st convolution: %.16f \n',err1);
    end
    Appp_test = conv(Ap,App_test,'same');
    dbf2 = figure("Name",'dbf1','Visible','off');
    plot(length(Ap)/2:length(Ap),Appp(length(Ap)/2:length(Ap))); hold on; 
    plot(length(Ap)/2:length(Ap),Appp_test(length(Ap)/2:length(Ap)),'--');
    err2 = abs(norm(Appp_test-Appp)/norm(Appp_test));
    if err2 > 10*eps %-> idem
       fprintf('Error on 2nd convolution: %.16f \n',err2);
    end
  end
end


