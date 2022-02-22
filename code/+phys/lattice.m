function gloc = lattice(zeta,D,lattice)
%% SELF-CONSISTENCY relation for the DMFT loop
%  
%  Most generally the DMFT self-consistency relation expresses the local
%  GF as a the Fourier transform of the lattice GF evaluated at r=0, namely
%  \sum_k G_latt and the self-energy appearing in this expression should be
%  substituted with the impurity one (locality assumption), thus giving:
%
%          Gloc(z) = F_lattice[Sloc(z)],             DOS <--> lattice
%
%                  = H[DOS](z-Sloc(z)),   i.e. F_lattice <--> H[DOS]
%
%  where H[u](:) is the HILBERT transform of the function u(..)
%
%  > On the BETHE lattice we can avoid computing numerically with hilbert()
%    for that the DOS is very simple (semi-circular of radius D) and the 
%    transform can be written explicitly as:
%
%            H(z) = 2(z-sign[Re(z)]sqrt[z^2-D^2])/D^2
%
%                 = 2(z-sign[Im(s)]*s)/D^2, with s = sqrt[z^2-D^2]
%
%    where the last passage is a trick to allow U=0 (Sloc(:)=0) evaluations 
%
%% Theoretical Background at:
%
%  http://www.physics.rutgers.edu/~haule/681/Perturbation.pdf
%
%  https://en.wikipedia.org/wiki/Hilbert_transform
%
%% BSD 3-Clause License
%
%  Copyright (c) 2020-2022, Gabriele Bellomia
%  All rights reserved.

    switch lattice
        
        case 'bethe'

            s    = sqrt(zeta.^2 - D^2);
            p    = sign(imag(s)) .* s;
            gloc = 2 * (zeta - p) / D^2;
            
        case 'square'
            
            invz = 1./zeta;
            ellk = ellipticK(D^2*invz.^2);
            gloc = 2/pi.*invz.*ellk;
            
        case 'chain'
            
            invz = 1./zeta;
            fact = D.*invz;
            gloc = fact./sqrt(1-fact.^2);
            
        otherwise
            
            error('Invalid lattice');
        
    end
    
end

%% ellipticK(m) evaluates to:
%
%       π
%       ─
%       2
%       ⌠
%       ⎮         1
%       ⎮ ────────────────── dϕ
%       ⎮    _______________
%       ⎮   ╱          2
%       ⎮ ╲╱  1 - m⋅sin (ϕ)
%       ⌡
%       0
%
% NB: it requires Symbolic Math Toolbox! For real m we could use the faster 
%     and built-in ellipke(m), but here we unfortunately need { m ∈ ℂ }
%
