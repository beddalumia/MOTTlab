function line = Tline(Z,U,Tmin,Tstep,Tmax)
    line = figure("Name",'T-driven MIT');
    T = Tmin:Tstep:Tmax;
    scatter(T,Z,'k','filled');
    title(sprintf('IPT  |  Quasiparticle weight at U/D = %f',U))
    xlabel('$T$','Interpreter','latex')
    ylabel('$Z = M/M^*$','Interpreter','latex')
end


