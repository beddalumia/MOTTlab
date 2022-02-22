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
            ellk = elliptic(D^2*invz.^2,lattice);
            gloc = 2/pi.*invz.*ellk;
            
        case {'cubic','sc'}
            
            invd = 3/D;
            zeta = invd.*zeta;
            zsqr = zeta.^(-2);
            xi   = sqrt(1-sqrt(1-zsqr)) ./ sqrt(1+sqrt(1-9.*zsqr));
            invx = 1 ./ ((1 - xi).^3 .* (1 + 3*xi));
            ellk = elliptic(16 .* xi.^3 .* invx, lattice);
            gloc = invd .* (1-9.*xi.^4) .* (2/pi.*ellk).^2 .* invx ./ zeta;
            
        case 'bcc'
            
            zren = zeta/D;
            ellk = elliptic(0.5 .* (1-sqrt(1-zren.^(-2))), lattice);
            gloc = 4 ./ (pi^2 .* zeta) .* ellk.^2; 
            
        case 'chain'
            
            invz = 1./zeta;
            fact = D.*invz;
            gloc = fact./sqrt(1-fact.^2);
            
        otherwise
            
            error('Invalid lattice');
        
    end
    
end

function k = elliptic(m,str)
%% elliptic(m,str) wraps ellipticK(m) which in turn evaluates:
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
fprintf(2,'SLOW-RUN: %s lattice requires symbolic computations \n',str);
k = ellipticK(m);
end
