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
%             Half-bandwidth of the noninteracting system
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
    f = phys.fermi(w,beta);
    
    counter = 0;
    CONVERGED = false;
    LOOP = true;
    
    g0 = gloc; % amounts to sloc=0 initizalization
    
    while LOOP
        
        % Increment the loop counter
           counter = counter + 1;
            
        % Spectral-function of Weiss field
            A0 = -imag(g0) ./ pi;
            
        % Enforce particle-hole and half-filling ( -> stability )
            A0 = 0.5 * (A0 + flip(A0)); % [ flip{v(1:end)}=v(end:1) ]
            
        % Second Order Perturbation Theory, to get Im(Sigma(w))
            imsloc = phys.sopt(A0,f,U)*dw^2;
            
        % Kramers-Kronig transform, to get Re(Sigma(w))
            sloc = math.fkkt(imsloc) + 1i.*imsloc;
            
        % Store "old" fields
            old_gloc = gloc; old_g0 = g0;
            
        % Self-Consistency relations
            gloc = phys.gloc(w-sloc,D,dos);  
            g0   = 1./(1./gloc + sloc);
            
        % Mixing ( -> stability )
            x  = old_g0;
            Fx = g0 - old_g0;
            g0 = adaptive_mixing(x,Fx,mix,counter);
            
        % Debug plotting
            %plot(w,real(g0)); hold on
            %plot(w,imag(g0)); hold off
            %pause
            
        % Logical Update 
            E = norm(gloc-old_gloc)/norm(old_gloc);
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

function x = adaptive_mixing(x,Fx,alpha,iter)
% ADAPTIVE MIXING SCHEME: directly imported from SciFortran
% > https://github.com/QcmPlab/SciFortran/tree/master/src/SF_OPTIMIZE

    alpha_max = 1;        % maximal mixing parameter α
    persistent Fx_prev    % stored Fx = (x_new - x)
    persistent beta       % stored β (adaptive mixing)

    assert(0<alpha && alpha<1, 'MIXING parameter must be in [0,1]')
    
    N = length(x);

    if(iter==1)
        Fx_prev = zeros(1,N);   % reset values at
        beta = alpha*ones(1,N); % first dmft loop
    end

    x = x + beta.*Fx; % = (1-β)x + βx_new

     if (iter > 1)
        % NATURAL LOOP VERSION
%         for j=1:N 
%             if(Fx_prev(j) * Fx(j) > 0) % If going smooth...
%                 beta(j) = beta(j) + alpha; % ...increase update speed
%                 if (beta(j) > alpha_max) 
%                     beta(j) = alpha_max; % but with a cutoff!
%                 end
%             else
%                 beta(j) = alpha;
%             end
%         end
        % VECTORIZED VERSION
        oscillating = real(Fx_prev).*real(Fx)<0 | imag(Fx_prev).*imag(Fx)<0;
        beta(not(oscillating)) = beta(not(oscillating)) + alpha;
        beta(beta>alpha_max) = alpha_max;
    end

    Fx_prev = Fx; % store persistent value

end
