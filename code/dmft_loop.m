function [gloc,sloc] = dmft_loop(gloc,w,U,beta,D,dos,mloop,mix,err,pmode)
%% DMFT_LOOP on a single-band Hubbard model at half filling. Using IPT.
%
%%    Parameters
%     ----------
%     gloc  : complex 1D ndarray
%             Local Green's function to use as seed/restart
%     w     : real 1D ndarray
%             Real frequency points
%     U     : float
%             Onsite interaction, Hubbard U
%     beta  : float
%             Inverse temperature
%     D     : float
%             Bandwidth of the noninteracting density of states
%     dos   : string [default: bethe]
%             Lattice on which the Hubbard model is defined
%     mloop : int
%             Max amount of DMFT iterations to perform
%     mix   : float \in [0,1]
%             Mixing parameter for iteration updates
%     err   : float \in [0,1]
%             Convergence threshold for self-consistency
%     pmode : string, optional [default: ~'quiet']
%             If 'quiet' deactivates loop-by-loop prints
% 
%%    Returns
%     -------
%     gloc  : complex 1D ndarray
%             DMFT iterated (converged) local Green's function
%     sloc  : complex 1D ndarray
%             DMFT iterated (converged) self-energy
%
%% Theoretical Background at:
%
%  http://www.physics.rutgers.edu/~haule/681/Perturbation.pdf
%
%% BSD 3-Clause License
%
%  Copyright (c) 2020, Gabriele Bellomia
%  All rights reserved.

if(~exist('pmode','var'))
    pmode = 'notquiet';
end
quiet = strcmp(pmode,'quiet');

    %% Iterated Perturbation Theory (IPT)

    dw = w(2)-w(1);
    eta = 2i * dw;
    f = phys.fermi(w,beta);
    counter = 0;
    CONVERGED = false;
    LOOP = true;

    while LOOP
        
        % Weiss field from local Green's function
            g0 = 1 ./ (w + eta - 0.25 .* gloc);
            
        % Spectral-function of Weiss field
            A0 = -imag(g0) ./ pi;
            
        % Enforce particle-hole and half-filling ( -> stability )
            A0 = 0.5 * (A0 + flip(A0)); % [ flip{v(1:end)}=v(end:1) ]
            
        % Second Order Perturbation Theory, to get Im(Sigma(w))
            imsloc = phys.sopt(A0,f,U)*dw^2;
            
        % Kramers-Kronig transform, to get Re(Sigma(w))
            sloc = math.fkkt(imsloc) + 1i.*imsloc;
            
        % Self-Consistency relation
            new_gloc = phys.lattice(w-sloc,D,dos); % D is the dos bandwidth
            
        % Mixing ( -> stability )
            old_gloc = gloc; gloc = mix*new_gloc+(1-mix)*old_gloc; 
            
        % Logical Update
            counter = counter + 1;
            E = norm(gloc-old_gloc)/norm(gloc);
            if E < err
               CONVERGED = true; 
            end
            LOOP = (counter < mloop) && not(CONVERGED);
            
        % Print Info
            if(~quiet)
            fprintf('DMFT-loop #%d (max %d) has ended\n', counter, mloop);
            end
            
    end
    
    if(~quiet)
        fprintf('\n');
    end
    
    if CONVERGED
        fprintf(1,'DMFT has converged after %d steps\n', counter);
    else
        beep
        fprintf(2,'DMFT *not* converged after %d steps\n', counter);
    end
    
    fprintf('> error = %f\n\n',E);
    
end




