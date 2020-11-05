%% BSD 3-Clause License
% 
% Copyright (c) 2020, Gabriele Bellomia
% All rights reserved.

function [gloc,Sigma] = DMFT_loop(gloc,w,D,U,beta,loops,mixing,error)
%% DMFT_LOOP Single band Bethe lattice at half filling. Using IPT.
%
%%    Parameters
%     ----------
%     gloc : complex 1D ndarray
%         local Green's function to use as seed
%     w : real 1D ndarray
%         real frequency points
%     D : float
%         Radius of the semicircular Density of States, Bandwidth
%     U : float
%         Onsite interaction, Hubbard U
%     beta : float
%         Inverse temperature
%     loops : int
%         Max amount of DMFT iterations to perform
%     mixing : float \in [0,1]
%         Mixing parameter for iteration updates
%     error : float \in [0,1]
%         Convergence threshold for self-consistency
% 
%%    Returns
%     -------
%     gloc : complex 1D ndarray
%         DMFT iterated (converged) local Green's function
%     Sigma : complex 1D ndarray
%         DMFT iterated (converged) self-energy

%% Iterated Perturbation Theory (IPT)

dw = w(2)-w(1);
eta = 2i * dw;
f = FermiDirac(w, beta);
counter = 0;
SelfCons = false;
ToLoop = true;

    while ToLoop == true
        % Self-consistency
            g0 = 1 ./ (w + eta - 0.25 .* gloc);
        % Spectral-function of Weiss field
            A0 = -imag(g0) ./ pi;
        % Clean for PH and Half-fill (????)
            A0 = 0.5 * (A0 + flip(A0)); %( flip([1 2 3]) == [3 2 1] )%
        % 2nd order Perturbation Theory (PT)
            isi = PT(A0,f,U)*dw^2;
        % Kramers-Kronig relation, using built-in FFT (hilbert subroutine)
            Nyquist = length(isi)*4;
            H = hilbert(isi,Nyquist);
            hsi = imag(-H(1:length(isi)));
            Sigma = hsi + 1i * isi;
        % Semicircular Hilbert Transform ( != hilbert subroutine -> self-coded )
            new_gloc = semiHilbert(w-Sigma,D);
        % Mixing (for convergence stability purposes)
            old_gloc = gloc;
            gloc = mixing*new_gloc+(1-mixing)*old_gloc; % D is the DOS "radius"
        % Logical Update
            counter = counter + 1;
            E = norm(gloc-old_gloc)/norm(gloc);
            if E < error
               SelfCons = true; 
            end
            ToLoop = (counter < loops) && (SelfCons == false);
    end
    if SelfCons == true
        fprintf('DMFT has converged after %d steps\n', counter);
    else
        fprintf('DMFT has *not* converged after %d steps\n', counter);
        fprintf('[error = %f]\n',E);
    end
%% Theoretical Background at:
%  http://www.physics.rutgers.edu/~haule/681/Perturbation.pdf
end

