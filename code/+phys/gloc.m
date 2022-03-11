function gloc = gloc(zeta,D,lattice)
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
%% References
%
%  [gftools]    https://gftools.readthedocs.io/en/stable/index.html
%  [delves]     https://doi.org/10.1006/aphy.2001.6148
%  [kogan]      https://doi.org/10.4236/graphene.2021.101001
%  [horiguchi]  https://doi.org/10.1063/1.1666155
%  [morita]     https://doi.org/10.1063/1.1665693
%
%% BSD 3-Clause License
%
%  Copyright (c) 2020-2022, Gabriele Bellomia
%  All rights reserved.

    switch lattice
        
        case 'bethe'

            gloc = gloc_bethe(zeta,D);
            
        case {'square','2d'}
            
            gloc = gloc_square(zeta,D);
            
        case {'cubic','3d','sc'}
            
            gloc = gloc_cubic(zeta,D);
            
        case 'bcc'
            
            gloc = gloc_bcc(zeta,D);
            
        case {'chain','1d'}
            
            gloc = gloc_chain(zeta,D);
            
        case 'lieb'
            
            gloc = gloc_lieb(zeta,D);
            
        case {'honey','honeycomb','graphene'}
            
            gloc = gloc_honey(zeta,D);
            
        otherwise
            
            error('Invalid lattice');
        
    end
    
end

function gloc = gloc_bethe(zeta,D)
% Bethe lattice (HM-∞d)

    s    = sqrt(zeta.^2 - D^2);
    p    = sign(imag(s)) .* s;
    gloc = 2 * (zeta - p) / D^2;

end

function gloc = gloc_chain(zeta,D)
% Chain (HM-1d)

    invz = 1./zeta;
    fact = D.*invz;
    gloc = invz./sqrt(1-fact.^2);
            
end

function gloc = gloc_square(zeta,D)
% Square lattice (HM-2d)
% Ref. to [kogan]

    invz = 1./zeta;
    ellk = elliptic(D^2*invz.^2);
    gloc = 2/pi.*invz.*ellk;

end

function gloc = gloc_cubic(zeta,D)
% Simple-Cubic lattice (HM-3d)
% Ref. to [delves]

    invd = 3/D;
    zeta = invd.*zeta;
    zsqr = zeta.^(-2);
    xi   = sqrt(1-sqrt(1-zsqr)) ./ sqrt(1+sqrt(1-9.*zsqr));
    invx = 1 ./ ((1 - xi).^3 .* (1 + 3*xi));
    ellk = elliptic(16 .* xi.^3 .* invx);
    gloc = invd .* (1-9.*xi.^4) .* (2/pi.*ellk).^2 .* invx ./ zeta;

end

function gloc = gloc_bcc(zeta,D)
% Body-Centered-Cubic lattice
% Ref. to [morita]

    zren = zeta/D;
    ellk = elliptic(0.5 .* (1-sqrt(1-zren.^(-2))));
    gloc = 4 ./ (pi^2 .* zeta) .* ellk.^2; 
            
end

function gloc = gloc_lieb(zeta,D)
% Lieb lattice (flat-band metal)
% Ref. to [kogan]

    dren = D/2^1.5;
    zren = zeta/dren;
    peak = 1/3./zeta;
    sqlt = gloc_square(zren.^2 - 4, 4);
    gloc = peak + 2/3*zren/dren .* sqlt;

end

function gloc = gloc_honey(zeta,D)
% Honeycomb (graphene) lattice
% Ref. to [horiguchi]

    gloc = gloc_hexa(zeta,D);

end

function gloc = gloc_hexa(zeta,D)
% Hexagonal (triangular) lattice
% Ref. to [horiguchi]

    dren = D*4/9;
    zren = zeta/dren;
    sing = (zren * dren == -1); zren(sing) = 0;
    rr   = csqrt(2*zren + 3);
    gg   = 4 ./ (csqrt(rr - 1).^3 .* csqrt(rr + 3));
    kk   = csqrt(rr) .* gg;
    mm   = kk.^2;
    ellk = elliptic(mm);
    ikp  = imag(kk) > 0;
    ellk(ikp) = ellk(ikp) + 2i * elliptic(1 - mm(ikp));
    gloc = 1/(pi*dren) .* gg .* ellk;
    gloc(sing) = -1i*inf;
    
    function s = csqrt(z)
    % appropriate branch of sqrt for the triangular lattice
        s = ones(size(z));
        s( real(z) < 0 & imag(z) < 0 ) = -1; % sign switch
        s = s .* sqrt(z);
    end

end

function k = elliptic(m)
% ELLIPTIC(m) wraps a homebrew variant of ellipke allowing a numerical*
% evaluation of the complete elliptic integral of the first kind:
%
%                 π
%                 ─
%                 2
%                 ⌠
%                 ⎮         dϕ
%                 ⎮ ────────────────── ,  >>> for m ∈ ℂ <<<
%                 ⎮    _______________
%                 ⎮   ╱          2
%                 ⎮ ╲╱  1 - m⋅sin (ϕ)
%                 ⌡
%                 0
%
% *It can be benchmarked with a reliable (but severely slower) symbolic 
%  implementation provided by MathWorks within the Symbolic Math Toolbox,
%  namely ellipticK(m). Here we remove any explicit test to keep the DEBUG 
%  mode actually runnable and avoid any needless dependency on toolboxes.
%
% See also: math.cellke, ellipke, ellipticK

  k = math.cellke(m);       
  
end
