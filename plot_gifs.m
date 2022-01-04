% Assumes main_Bethe_dmft has been run (ugly but for now that's it)

if Uline
    %% U-driven MIT span: spectral gifs
    U = 0; i = 0;
    while U <= 5.00 
        i = i + 1;
        Sigma = Sigma_loc{i};
        G = gloc{i};
        Z = Zweight(w,Sigma);
          % Uncomment for gifs of the span
            TitleString = sprintf('U=%f', U)
            f = figure("Name",'Spectral Function (DOS)')
            FilledStates = w(w<=0);
            FilledDOS = -imag(G(w<=0));
            EmptyStates = w(w>0);
            EmptyDOS = -imag(G(w>0)); 
            area(FilledStates, FilledDOS, 'FaceColor', [0.7 0.7 0.7]); hold on
            area(EmptyStates, EmptyDOS, 'FaceColor', [1 1 1]);
            title(sprintf('IPT  |  DOS @ U/t = %.2f, beta = %d', 2*U, beta));
            xlabel('\omega')
            ylabel('\pi A(\omega)')
            ylim([0,2]);
            legend('Filled States','Empty States')
            print(append('DOS\',TitleString,'.png'),'-dpng','-r1200')
            close(f);
            f = figure("Name",'Self-Energy')
            plot(w,real(Sigma)./U,'LineWidth',1); hold on
            plot(w,imag(Sigma)./U,'LineWidth',1);
            title(sprintf('IPT  |  Self-Energy @ U/t = %.2f, beta = %d', 2*U, beta));
            xlabel('\omega')
            ylabel('Units of U');
            ylim([-1,1]);
            legend('Re \Sigma_{loc}(\omega)','Im \Sigma_{loc}(\omega)')
            print(append('Sigma\',TitleString,'.png'),'-dpng','-r1200')
            close(f);
        U = U + Ustep;
    end 
end 



    