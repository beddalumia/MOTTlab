function line = Uline(Z,I,beta,Umin,Ustep,Umax,D)
    line = figure("Name",'U-driven MIT');
    U = Umin:Ustep:Umax;
    scatter(U/D,Z,'k','filled'); hold on
    scatter(U/D,I,'r','filled');
    title(sprintf('IPT  |  U-driven MIT at T = %f',1/beta))
    xlabel('$U/D$','Interpreter','latex')
    legend({'$Z_F = 1 / ( 1 - Re[\partial_\omega\Sigma(0)] )$',...
        '$I_L = 1/\pi\int d\omega Im[G(\omega)\partial_\omega\Sigma(\omega)]$'},...
        'Interpreter','latex')
end


