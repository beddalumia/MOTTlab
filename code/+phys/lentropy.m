function lS = lentropy(n,D)
    
    % Local entropy (Mod.Phys.Lett.B.27:05)
    lS = -2*(n/2-D)*log2(n/2-D)-2*D*log2(D);

end