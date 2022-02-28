function line = Tline(Z,U,Tmin,Tstep,Tmax,D)
    line = figure("Name",'T-driven MIT');
    T = Tmin:Tstep:Tmax;
    scatter(T,Z,'k','filled'); box on;
    title(sprintf('IPT  |  Quasiparticle weight at U/D = %f',U/D))
    xlabel('$T$','Interpreter','latex')
    ylabel('$Z = M/M^*$','Interpreter','latex')
end


