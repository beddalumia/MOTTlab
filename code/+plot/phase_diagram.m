function map = phase_diagram(marker,Umin,Ustep,Umax,Tmin,Tstep,Tmax,D)
    X = Umin:Ustep:Umax; X=X/D;
    Y = Tmin:Tstep:Tmax;
    Z = marker;
    % Surface plot
    map = figure("Name",'Phase Diagram');
    surf(X,Y,Z); colormap('winter')
    title(sprintf('IPT  |  (U,T) Phase-Diagram'))
    xlabel('$U/D$','Interpreter','latex')
    ylabel('$T$','Interpreter','latex')
end

