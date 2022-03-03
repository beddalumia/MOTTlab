function docc = docc(iw,sloc,gloc,dens,u)

    docc = 0.25;
    if u > 0
       beta = pi/iw(1);
       epot = 2/beta * sum(real(sloc.*gloc));
       ehar = (0.5-dens)*u/2;
       docc = (epot-ehar)/u;
    end

end