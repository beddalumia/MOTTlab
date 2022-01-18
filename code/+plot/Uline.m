function line = Uline(Z,beta,Umin,Ustep,Umax)
    line = figure("Name",'U-driven MIT');
    U = Umin:Ustep:Umax;
    scatter(2.*U,Z,'k','filled'); hold on
    title(sprintf('IPT  |  U-driven MIT at T = %f',1/beta))
    xlabel('$U/t$','Interpreter','latex')
    ylabel('$Z = M/M^*$','Interpreter','latex')
end


