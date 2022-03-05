function docc = docc(iv,smatsu,gmatsu,uloc)
%% Computes double occupancy ⟨c^†_up c_up c^†_dw c_dw⟩ via Matsubara sums
%
%  Input:
%       iv      : real valued array, fermionic frequencies, ν = (2n+1)πT
%       smatsu  : complex valued array, matsubara self-energy, ∑(iν)
%       gmatsu  : complex valued array, matsubara green's function, G(iν)
%       uloc    : real double, local (on-site) Hubbard interaction, U
%   Output:
%       docc    : real double, double occupancy, D = ⟨n_up n_dw⟩
%
%   The implementation is based on Matsubara formalism, as briefly reported
%   in [PRB.93.155162] (which shows how to go beyond this simple approach). 
%
%   The idea is to express the potential energy as 
%
%                         Epot = ⟨H_int⟩ = U⟨n_up n_dw⟩
%
%   so that the double occupancy is just given by Epot/U. In turn we can
%   evaluate the potential energy by standard finite temperature formalism,
%   as given by [Fetter-Walecka] in eq. 23.14, by trasforming to imaginary
%   frequency and assuming ∑(k,iν) to be local (which amounts to DMFT).
%
%   This would simply give:
%
%                         Epot = 2/beta \sum_iν ∑(iν)G(iν)
%
%   But needs to be corrected by the Hartree term (the formula assumes an
%   Hubbard interaction of the form H_int ~ U n_up n_dw, but our convention
%   that half-filling corresponds to µ = 0 requires instead to have:
%
%                 ⟨H_int⟩ = ⟨U(n_up-1/2)(n_dw-1/2)⟩
%                         = U * ⟨n_up n_dw⟩ - ⟨n_up/2 + n_dw/2 - 1/4⟩ * U 
%                         = U * D - (n - 1/2) * U/2 
%
%                    ⟹ D = (Epot - U(1/2 - n)/2) / U 
%
%   The noninteracting (U = 0) limit has to be treated separately, for that
%   we cannot divide by zero. The result is trivial since for free fermions
%   the two spin polarizations are totally independent and we have that:
%
%             ⟨n_up n_dw⟩ = ⟨1 1⟩ = ⟨1 0⟩ = ⟨0 1⟩ = ⟨0 0⟩ = 0.25
%
%   References:
%
%       [PRB.93.155162]  = Double occupancy in dynamical mean-field theory 
%                          and the dual boson approach, E. van Loon et al.
%       [Fetter-Walecka] = Quantum Many-Particle Systems, Dover 2003
%
%% BSD 3-Clause License
%
%  Copyright (c) 2022, Gabriele Bellomia
%  All rights reserved.

    % Density (hardcoded to half-filling)
    dens = 1.00;
    
    % Noninteracting double occupancy
    docc = 0.25;
    
    % Interacting double occupancy
    if uloc > 0
       beta = pi/iv(1);
       epot = 2/beta * sum(real(smatsu.*gmatsu));
       ehar = (0.5 - dens)/2 * uloc;
       docc = (epot - ehar)  / uloc;
    end
    
end

