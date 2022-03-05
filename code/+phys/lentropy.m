function S = lentropy(iv,smatsu,gmatsu,uloc,dvec)
%% Computes the local entropy, via Matsubara summation
%
%  Input:
%       iv      : real valued array, fermionic frequencies, ν = (2n+1)πT
%       smatsu  : complex valued array, matsubara self-energy, ∑(iν)
%       gmatsu  : complex valued array, matsubara green's function, G(iν)
%       uloc    : real double, local (on-site) Hubbard interaction, U
%       dvec    : real valued array, OPTIONALLY provided double occupancy
%   Output:
%       sloc    : real double, local entropy, S = -Tr(ρ_loc*log2(ρ_loc)) 
%
%  The local density matrix ρ_loc can be defined directly considering the
%  known 4 pure states of a single-site many-fermion system (the Anderson
%  impurity in the DMFT scheme).
%
%            ⎛⟨(1-n_up)(1-n_dw)⟩      0             0            0     ⎞
%    ρ_loc = ⎜        0         ⟨n_up(1-n_dw)⟩      0            0     ⎟
%            ⎜        0               0       ⟨(1-n_up)n_dw⟩     0     ⎟
%            ⎝        0               0             0       ⟨n_up n_dw⟩⎠
%
%  Hence:
%
%        S = -(1 - n_up - n_dw + D) * log2(1 - n_up - n_dw + D)
%
%          = -2 * (n/2 - D) * log2(n/2 - D)  -  2 * D * log2(D)
% 
%  where n = ⟨n_up + n_dw⟩ and D = ⟨n_up n_dw⟩ are respectively the local
%  density and the double occupancy. The former is hardcoded to 1.d0 here,
%  for we are enforcing half-filling everywhere in the codebase, while the
%  latter is retrieved via a suitable Matsubara sum, by calling phys.docc()
%
%  To avoid wasting time performing the Matsubara sums twice, you can feed
%  phys.entropy() with an optional vector of already computed D values, and
%  so skip the call to phys.docc, for a fast vectorized evaluation of S. :)
%
%  See also phys.docc
%
%% BSD 3-Clause License
%
%  Copyright (c) 2022, Gabriele Bellomia
%  All rights reserved.
    
    % Density (hardcoded to half-filling)
    n = 1.00;

    % Double occupancy (Phys.Rev.B.93.155162)
    if nargin < 5
       D = phys.docc(iv,smatsu,gmatsu,uloc);
    else
       D = dvec;   
    end
    
    % Local entropy (Mod.Phys.Lett.B.27:05)
    S = -2*(n/2-D).*log2(n/2-D)-2*D.*log2(D);

end



