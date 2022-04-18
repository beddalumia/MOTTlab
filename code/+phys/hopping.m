function t = hopping(D,LATTICE)
                                                              global DEBUG

    switch LATTICE

        case {'bethe','chain','1d'}

            t = D/2;

        case {'2d-toy','rect'}

            t = NaN;

if DEBUG
            warning("We want to emulate a generic " + ...
            "2d material, but the DOS does not really " + ...
            "correspond to any lattice model, so the t(D) " + ...
            "relation is totally arbitrary. To be investigated. ")
end

        case {'square','2d'}

            t = D/4;

        case {'sc','cubic','3d'}

            t = D/6;

        case 'bcc'

            t = D/8;

        case {'honey','honeycomb','graphene'}

            t = D/3;
            
        case 'lieb'
            
            t = D*2^1.5;

        otherwise

            error("Invalid lattice: " + ...
            "see 'help phys.gloc' for the available choices.");

    end

end