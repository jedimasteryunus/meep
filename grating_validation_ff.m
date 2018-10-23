function grating_validation_ff()
    fname = 'grating_validation.out';
    
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
    
    ang = asin(55/800);
    
    plot(angs, Px, angs, Py, angs, Pz)
    
    figure

    if true
        polar(angs, Pr/max(Pr))
        hold on

        center = pi/2 + asin(55/800);

        polar([center + ang, 0, center - ang], [1 0 1])
        hold off

        sum(Pr(angs <= center + ang & angs >= center - ang))/sum(Pr)
        sum(Pr(~(angs <= center + ang & angs >= center - ang)))/sum(Pr)
    else
        plot(angs, Pr/max(Pr))
    end
end