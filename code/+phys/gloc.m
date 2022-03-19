function gloc = gloc(zeta,D,lattice)
%% Local GF for a lattice with D half-bandwidth, evaluated at z ∈ ℂ
%
%  Input:
%       zeta    : complex valued array, domain on which to evaluate the GF
%       D       : real double, hald-bandwidth of the noninteracting system
%       lattice : string, to select the underlying lattice [default: bethe]
%  Output:
%       gloc    : complex valued array, local GF, G(z) = ∫dε DOS(ε)/(z-ε)
%
%  You may recognize the integral expression of G(z) as the Hilbert transform
%  of the noninteracting spectral function A₀(ω), aka density of states DOS,
%  evaluated at complex argument z. (see [hilbert])
%  But for proper use within a DMFT code it is *crucial* that z could really
%  take *any* complex value, not just those on real- or imaginary-axis. This
%  makes impossible to make use of the efficient hilbert() function provided
%  within the Signal Processing Toolbox, as for the Kramers Kronig transform
%  implemented in math.fkkt(). The reason for this very general requirement
%  lies in the main application within DMFT: the self-consistency relation.
%
%  Most generally the DMFT self-consistency relation expresses the local GF
%  as the lattice GF evaluated at r = 0, namely ∑ₖG(k,ω) = ∑ₖ1/(ω-εₖ-∑(k,ω))
%  where the self-energy appearing should be assumed k-independent, hence 
%  substituted with the impurity one (locality assumption), thus giving:
%
%      ∑ₖG(k,ω) = ∑ₖ1/(ω-εₖ-∑(ω)) = ∫dε DOS(ε)/(ω-ε-∑(ω)) = G(ω-∑(ω))
%
%               = H[DOS](ω-∑(ω)),
%
%  where H[u](z) is the HILBERT transform of the function u(:), evaluated
%  at the most generic complex argument z = ω-∑(ω). See also [haule_ln].
%
%  To secure performance we *must* avoid computing such hilbert transform
%  by direct integration, either by sum/trapz or adaptive quadrature, thus
%  exploiting instead lattice-specific (semi)analytical insight.
%
%  The default choice for DOS(ε) corresponds to the infinite-coordination
%  Bethe lattice ['bethe'], but many others model (semi)metals have been
%  added. You can request a lattice by passing one of the supported strings.
%
%  • On the BETHE lattice we can virtually avoid any numerical computation
%    for that the DOS is very simple (semicircle of radius D) and the
%    transform can be integrated explicitly as:
%
%           G(z) = 2(z-sign[Re(z)]sqrt[z^2-D^2])/D^2
%
%                = 2(z-sign[Im(s)]*s)/D^2, with s = sqrt[z^2-D^2]
%
%    where the last passage is a trick to allow U=0 (∑(:)=0) evaluations.
%
%    >> phys.gloc(zeta,D,'bethe')      ⟹ infinite dimensions
%
%  • The 1d CHAIN can be also evaluated by full analytical means, giving
%
%                   +π
%                   ⌠
%                   ⎮  dϕ          1                    1         
%           G(z) =  ⎮ ──── ⋅ ──────────────  =  ─────────────────
%                   ⎮  2π     z - D⋅cos(ϕ)      z⋅sqrt(1-(D/z)^2)
%                   ⌡
%                  -π
%
%    as proved in [economou].
%
%    >> phys.gloc(zeta,D,'chain')      or
%    >> phys.gloc(zeta,D,'1d')         ⟹ 1d chain
%
%  • The 2d and 3d HYPERCUBIC lattices can be recasted as suitable elliptic 
%    integrals(†), for more info see [economou], [delves], [morita].
%
%    >> phys.gloc(zeta,D,'square')     or
%    >> phys.gloc(zeta,D,'2d')         ⟹ 2d square lattice
%    >> phys.gloc(zeta,D,'cubic')      or
%    >> phys.gloc(zeta,D,'sc')         or
%    >> phys.gloc(zeta,D,'3d')         ⟹ 3d simple-cubic lattice
%    >> phys.gloc(zeta,D,'bcc')        ⟹ 3d body-centered-cubic
%
%  • The 2d LIEB lattice can be expressed in terms of a 2d-square lattice and
%    a naive ~1/z pole representing the flat band at the Fermi level. [kogan]
%
%    >> phys.gloc(zeta,D,'lieb')       ⟹ 2d lieb lattice model
%
%  • The 2d HONEYCOMB lattice can be expressed in terms of a triangular cell,
%    which in turn evaluates _two_ complete elliptic integrals(†). Hence this
%    is the numerically heaviest local Green's function implemented so far.
%
%    See [horiguchi] for all details.
%
%    >> phys.gloc(zeta,D,'honeycomb')  or
%    >> phys.gloc(zeta,D,'honey')      or
%    >> phys.gloc(zeta,D,'graphene')   ⟹ 2d honeycomb structure
%
%  (†): The elliptic integrals are the main challenge here, for the standard
%       numerical implementations support only very restricted input domains
%       and so are not suited to the evaluation of G(z) at complex z values.
%       Here we provide a homebrew extension of the standard ellipke(), that
%       has been benchmarked against the symbolic ellipticK() provided by the
%       Symbolic Math Toolbox (which is reliable but of course very slow).
%
%  See also math.zellke
%
%% References:
%  [haule_ln]   http://www.physics.rutgers.edu/~haule/681/Perturbation.pdf
%  [hilbert]    https://en.wikipedia.org/wiki/Hilbert_transform
%  [gftools]    https://gftools.readthedocs.io/en/stable/index.html
%  [economou]   https://doi.org/10.1007/3-540-28841-4_5
%  [delves]     https://doi.org/10.1006/aphy.2001.6148
%  [kogan]      https://doi.org/10.4236/graphene.2021.101001
%  [horiguchi]  https://doi.org/10.1063/1.1666155
%  [morita]     https://doi.org/10.1063/1.1665693
%
%% BSD 3-Clause License
%
%  Copyright (c) 2020-2022, Gabriele Bellomia
%  All rights reserved.

    if nargin < 3

       lattice = 'bethe';

    end

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
            
            error("Invalid lattice: " + ...
            "see 'help phys.gloc' for the available choices.");
        
    end
    
end

function gloc = gloc_bethe(zeta,D)
% Bethe lattice (HM-∞d)
% Ref. to [economou]

    s    = sqrt(zeta.^2 - D^2);
    p    = sign(imag(s)) .* s;
    gloc = 2 * (zeta - p) / D^2;

end

function gloc = gloc_chain(zeta,D)
% Chain (HM-1d)
% Ref. to [economou]

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
% Ref. to [delves], eqs(1.24 — 1.26)

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
% Ref. to [morita], eqs(2.1, 2.4)

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
% Honeycomb (graphene) structure
% Ref. to [horiguchi]

    dren = D/1.5;
    zren = zeta/dren;
    hexa = gloc_hexa(2 * zren.^2 - 1.5, 9/4);
    gloc = 2/dren * zren .* hexa;

end

function gloc = gloc_hexa(zeta,D)
% Hexagonal (triangular) lattice
% Ref. to [horiguchi]

    dren = D*4/9;
    zren = zeta/dren;
    % † must exploit hermiticity
    zadv = (imag(zren) < 0);
    % † > compute only the retarded part 
    zren(zadv) = conj(zren(zadv));
    % ‡ must handle singularity: safe zero value
    sing = (zren * dren == -1); zren(sing) = 0;
    rvec = csqrt(2*zren + 3);
    % eq(2.9)
    gvec = 4 ./ (csqrt(rvec - 1).^3 .* csqrt(rvec + 3));
    % eq(2.11)
    kvec = csqrt(rvec) .* gvec;
    mvec = kvec.^2;
    ellk = elliptic(mvec);
    % eqs(2.22, 2.18)
    ikp  = imag(kvec) > 0; % fix correct branch
    ellk(ikp) = ellk(ikp) + 2i * elliptic(1 - mvec(ikp));
    % eq(2.6)
    gloc = 1/(pi*dren) .* gvec .* ellk;
    % ‡ proper singular dispatch: (-i∞)
    gloc(sing) = -1i*inf;
    % † > build the advanced by symmetry
    gloc(zadv) = conj(gloc(zadv));
    
    function s = csqrt(z)
    % Appropriate branch of sqrt for the triangular lattice
    % Ref. to [gftools]
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
%  namely ellipticK(m). Here we avoid any explicit test to keep the DEBUG 
%  mode actually runnable and remove any needless dependency on toolboxes.
%
% See also: math.zellke, ellipke, ellipticK

  k = math.zellke(m);       
  
end
