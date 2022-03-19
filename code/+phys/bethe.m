function gloc = bethe(zeta,D)
%% Local GF for the Bethe lattice with D half-bandwidth, evaluated at z ∈ ℂ
%
%  Input:
%       zeta    : complex valued array, domain on which to evaluate the GF
%       D       : real double, hald-bandwidth of the noninteracting system
%  Output:
%       gloc    : complex valued array, local GF, G(z) = ∫dε DOS(ε)/(z-ε)
%
%  You may recognize the integral expression of G(z) as the Hilbert transform
%  of the noninteracting spectral function A₀(ω), aka density of states DOS,
%  evaluated at complex argument z.
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
%  at the most generic complex argument z = ω-∑(ω).
%
%  To secure performance we *must* avoid computing such hilbert transform
%  by direct integration, either by sum/trapz or adaptive quadrature, thus
%  exploiting instead lattice-specific analytical insight.
%
%  > On the BETHE lattice we can virtually avoid any numerical computation
%    for that the DOS is very simple (semicircle of radius D) and the
%    transform can be integrated explicitly as:
%
%           G(z) = 2(z-sign[Re(z)]sqrt[z^2-D^2])/D^2
%
%                = 2(z-sign[Im(s)]*s)/D^2, with s = sqrt[z^2-D^2]
%
%    where the last passage is a trick to allow U=0 (∑(:)=0) evaluations. 
%
%% Theoretical Background at:
%
%  http://www.physics.rutgers.edu/~haule/681/Perturbation.pdf
%
%  https://en.wikipedia.org/wiki/Hilbert_transform
%
%% BSD 3-Clause License
%
%  Copyright (c) 2020, Gabriele Bellomia
%  All rights reserved.

    s    = sqrt(zeta.^2 - D^2);
    p    = sign(imag(s)) .* s;
    gloc = 2 * (zeta - p) / D^2;
    
end



