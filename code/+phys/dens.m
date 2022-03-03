function dens = dens(w,gloc)

    dos  = -imag(gloc)/pi;
    dw   = w(end)-w(end-1);
    dens = sum(dos*dw);

end