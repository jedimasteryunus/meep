function [DL, angles, powers] = gratingSweep()
    DL = -10:5:10;
    angles = DL*0;
    powers = DL*0;
    
    cores = 4;

    for ii = 1:cores:length(DL)
        jj = ii;
        
        str = 'source activate mp_test & ';
        
        while jj < length(DL) && jj < ii + 4
%             str = [str 'source activate mp_test; python grating_validation.py ' num2str(DL(ii)) ' & '];
            str = [str 'source grating_validation.sh ' num2str(DL(ii)) ' & '];
            jj = jj + 1;
        end    
        
        str = [str 'wait'];
        
        str
        
        system(str);
        
        jj = ii;
        
        while jj < length(DL) && jj < ii + 4
            [angles(jj), powers(jj)] = grating_validation_ff(DL(jj));
            jj = jj + 1;
        end
    end
    
    plot(DL, angles);
    yyaxis right
    plot(DL, powers);
end