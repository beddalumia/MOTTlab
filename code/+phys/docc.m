function docc = docc(iv,smatsu,gmatsu,uloc)

    % Density (hardcoded to half-filling)
    dens = 1.00;
    
    % Noninteracting double occupancy
    docc = 0.25;
    
    % Interacting double occupancy
    % > see Phys. Rev. B 93 155162
    if uloc > 0
       beta = pi/iv(1);
       epot = 2/beta * sum(real(smatsu.*gmatsu));
       ehar = (0.5 - dens)/2 * uloc;
       docc = (epot - ehar)  / uloc;
    end
    
end