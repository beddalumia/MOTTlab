function line = Tline(Z,U,Tmin,Tstep,Tmax)
    line = figure("Name",'T-driven MIT');
    T = Tmin:Tstep:Tmax;
    scatter(T,Z,'k','filled'); box on;
    title(sprintf('IPT  |  Quasiparticle weight at U/t = %f',2*U)) % (Units: D=2t)
    xlabel('$T$','Interpreter','latex')
    ylabel('$Z = M/M^*$','Interpreter','latex')
end


