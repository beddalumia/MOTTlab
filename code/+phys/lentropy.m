function S = lentropy(iv,smatsu,gmatsu,uloc)
    
    % Density (hardcoded to half-filling)
    n = 1.00;

    % Double occupancy (Phys.Rev.B.93.155162) 
    D = phys.docc(iv,smatsu,gmatsu,uloc);
    
    % Local entropy (Mod.Phys.Lett.B.27:05)
    S = -2*(n/2-D)*log2(n/2-D)-2*D*log2(D);

end