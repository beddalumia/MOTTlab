% Assumes main_Bethe_dmft has been run (ugly but for now that's it)

if PhaseDiagram == true && ZeroTempSpan == true
    %% Zero Temperature MIT span
    U = 0; i = 0;
    while U <= 5.00 
        i = i + 1;
        Z = Zweight(w,Sigma_loc{i});
        scatter(2*U,Z,'k','filled'); hold on % (Units: D=2t)
        U = U + Ustep;
    end
    title(sprintf('IPT  |  Quasiparticle weight at T = %f',1/beta))
    xlabel('$U/t$','Interpreter','latex')
    ylabel('$Z = M/M^*$','Interpreter','latex')
    
elseif PhaseDiagram == true && ZeroTempSpan == false
    %% Full Phase Diagram span
    i = 0; T = 10^(-6); 
    while T < 0.05
        i = i + 1;
        j = 0; 
        U = 0; 
        while U <= 5.00  
            j = j + 1; beta = 1/T;
            Strength(j)= U; Temperature(i) = 1/beta;
            Z(i,j) = Zweight(w,Sigma_loc{i,j});
            U = U + Ustep;
        end
        T = T + Tstep; 
    end
    % Remove machine zeros
    %Z(Z<10^(-2)) = 0;
    %figure("Name","Z Contour")
    surf(Strength,Temperature,Z);
    
else
    %% Single (U, beta) pair calculation
    figure("Name",'Spectral Function (DOS)')
    FilledStates = w(w<=0);
    FilledDOS = -imag(gloc(w<=0));
    EmptyStates = w(w>0);
    EmptyDOS = -imag(gloc(w>0)); 
    area(FilledStates, FilledDOS, 'FaceColor', [0.7 0.7 0.7]); hold on
    area(EmptyStates, EmptyDOS, 'FaceColor', [1 1 1]);
    title(sprintf('IPT  |  DOS @ U = %.2f, beta = %d', U, beta));
    xlabel('\omega')
    ylabel('\pi A(\omega)')
    legend('Filled States','Empty States')
    figure("Name",'Self-Energy')
    plot(w,real(Sigma_loc),'LineWidth',1); hold on
    plot(w,imag(Sigma_loc),'LineWidth',1);
    title(sprintf('IPT  |  Self-Energy @ U = %.2f, beta = %d', U, beta));
    xlabel('\omega')
    legend('Re \Sigma_{loc}(\omega)','Im \Sigma_{loc}(\omega)')
    Z = Zweight(w,Sigma_loc)
    
end