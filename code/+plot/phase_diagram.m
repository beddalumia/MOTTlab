function map = phase_diagram(marker,Umin,Ustep,Umax,Tmin,Tstep,Tmax)
    X = 2.*(Umin:Ustep:Umax); % (Units: D=2t)
    Y = Tmin:Tstep:Tmax;
    Z = marker;
    % Surface plot
    map = figure("Name",'Phase Diagram');
    surf(X,Y,Z); colormap('winter')
    title(sprintf('IPT  |  (U,T) Phase-Diagram'))
    xlabel('$U/t$','Interpreter','latex')
    ylabel('$T$','Interpreter','latex')
end

