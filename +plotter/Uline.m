function line = Uline(Z,beta,Umin,Ustep,Umax)
    line = figure("Name",'U-driven MIT')
    i=0; U = Umin;
    while U <= Umax 
        i = i + 1;
        scatter(2*U,Z(i),'k','filled'); hold on % (Units: D=2t)
        U = U + Ustep;
    end
    title(sprintf('IPT  |  Quasiparticle weight at T = %f',1/beta))
    xlabel('$U/t$','Interpreter','latex')
    ylabel('$Z = M/M^*$','Interpreter','latex')
end


