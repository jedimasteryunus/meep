function [center, power] = grating_validation_ff(dw, dl)
%     fname = 'grating_validation.out';
    fname = ['grating_validation-w=' num2str(dw) 'nm-dl=' num2str(dl) 'nm.out'];
    
    mydata = dlmread(fname, ',');
    
    Ex=mydata(:,2); Ey=mydata(:,3); Ez=mydata(:,4); 
    Hx=mydata(:,5); Hy=mydata(:,6); Hz=mydata(:,7);

    Ex=conj(Ex); Ey=conj(Ey); Ez=conj(Ez);
    
    Px=real((Ey .* Hz)-(Ez .* Hy)); 
    Py=real((Ez .* Hx)-(Ex .* Hz));
    Pz=real((Ex .* Hy)-(Ey .* Hx));
    
    Pr=sqrt((Py.^2)+(Pz.^2));
    
%     mydata(:,2)
    
    angs = real(mydata(:,2));
    
    ang = asin(50/800);

%     center = pi/2 + asin(55/800);

    power = 0;
    centerbase = pi/2;
    center = centerbase;
    
    for testcenter = 0:.01:pi
        curpower = sum(Pr(angs <= testcenter + ang & angs >= testcenter - ang))/sum(Pr);
        
        if curpower > power
            power = curpower;
            center = testcenter;
        end
    end
    
%     plot(angs, Px, angs, Py, angs, Pz)
    
    figure

%     if true
        subplot(1,2,1)
        polar(angs, Pr/max(Pr))
        hold on

        polar([center + ang, 0, center - ang], [1 0 1])
        polar([centerbase + ang, 0, centerbase - ang], [1 0 1])
        hold off
%         sum(Pr(~(angs <= center + ang & angs >= center - ang)))/sum(Pr)
%     else
        subplot(1,2,2)
        plot(angs, Pr/max(Pr))
        hold on

        plot([center + ang, center + ang, center - ang, center - ang], [1 0 0 1])
        plot([centerbase + ang, centerbase + ang, centerbase - ang, centerbase - ang], [1 0 0 1])
        hold off
        
        ylim([0, 1])
%     end
    
    center = (center - centerbase)*180/pi;
end




