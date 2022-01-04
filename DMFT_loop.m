%% BSD 3-Clause License
% 
% Copyright (c) 2020, Gabriele Bellomia
% All rights reserved.

function [gloc,Sigma] = DMFT_loop(gloc,w,D,U,beta,mloop,mix,err,quiet)
%% DMFT_LOOP Single band Bethe lattice at half filling. Using IPT.
%
%%    Parameters
%     ----------
%     gloc  : complex 1D ndarray
%             Local Green's function to use as seed/restart
%     w     : real 1D ndarray
%             Real frequency points
%     D     : float
%             Radius of the semicircular Density of States, Bandwidth
%     U     : float
%             Onsite interaction, Hubbard U
%     beta  : float
%             Inverse temperature
%     mloop : int
%             Max amount of DMFT iterations to perform
%     mix   : float \in [0,1]
%             Mixing parameter for iteration updates
%     err   : float \in [0,1]
%             Convergence threshold for self-consistency
%     quiet : logical
%             If true deactivates loop-by-loop prints
% 
%%    Returns
%     -------
%     gloc  : complex 1D ndarray
%             DMFT iterated (converged) local Green's function
%     Sigma : complex 1D ndarray
%             DMFT iterated (converged) self-energy

%% Iterated Perturbation Theory (IPT)

dw = w(2)-w(1);
eta = 2i * dw;
f = FermiDirac(w, beta);
counter = 0;
SelfCons = false;
DoLOOP = true;

    while DoLOOP == true
        % Self-consistency
            g0 = 1 ./ (w + eta - 0.25 .* gloc);
        % Spectral-function of Weiss field
            A0 = -imag(g0) ./ pi;
        % Enforce particle-hole and half-filling
            A0 = 0.5 * (A0 + flip(A0)); % flip([1 2 3]) == [3 2 1] 
        % 2nd order Perturbation Theory (SOPT)
            isi = SOPT(A0,f,U)*dw^2;
        % Kramers-Kronig relation, using built-in FFT (hilbert subroutine)
            Nyquist = length(isi)*4;
            H = hilbert(isi,Nyquist);
            hsi = imag(-H(1:length(isi)));
            Sigma = hsi + 1i * isi;
        % Semicircular Hilbert Transform ( != hilbert internal function )
            new_gloc = BetheHilbert(w-Sigma,D);
        % Mixing (for convergence stability purposes)
            old_gloc = gloc;
            gloc = mix*new_gloc+(1-mix)*old_gloc; % D is the DOS "radius"
        % Logical Update
            counter = counter + 1;
            E = norm(gloc-old_gloc)/norm(gloc);
            if E < err
               SelfCons = true; 
            end
            DoLOOP = (counter < mloop) && (SelfCons == false);
        % Print Info
            if(~quiet)
            fprintf('DMFT-loop #%d (max %d) has ended\n', counter, mloop);
            end

    end
    if SelfCons == true
        fprintf('\nDMFT has converged after %d steps\n', counter);
    else
        fprintf('\nDMFT has *not* converged after %d steps\n', counter);
    end
    fprintf('> error = %f\n',E);
%% Theoretical Background at:
%  http://www.physics.rutgers.edu/~haule/681/Perturbation.pdf
end
