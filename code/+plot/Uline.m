function line = Uline(Z,I,beta,Umin,Ustep,Umax)
    line = figure("Name",'U-driven MIT');
    U = Umin:Ustep:Umax;
    scatter(2.*U,Z,'k','filled'); hold on; box on;
    scatter(2.*U,I,'r','filled'); % (Units: D=2t)
    title(sprintf('IPT  |  U-driven MIT at T = %f',1/beta))
    xlabel('$U/t$','Interpreter','latex')
    legend({'$Z_F = 1 / ( 1 - Re[\partial_\omega\Sigma(0)] )$',...
        '$I_L = 1/\pi\int d\omega Im[G(\omega)\partial_\omega\Sigma(\omega)]$'},...
        'Interpreter','latex')
end


