function line = Tline(Z,U,Tmin,Tstep,Tmax)
    line = figure("Name",'T-driven MIT');
    i=0; T = Tmin;
    while T <= Tmax
        i = i + 1;
        scatter(T,Z(i),'k','filled'); hold on
        T = T + Tstep;
    end
    title(sprintf('IPT  |  Quasiparticle weight at U/t = %f',2*U)) % (Units: D=2t)
    xlabel('$T$','Interpreter','latex')
    ylabel('$Z = M/M^*$','Interpreter','latex')
end


