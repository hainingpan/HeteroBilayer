function re=DeltaTT(h1,h2,parameters)
    %Delta_T^dagger in Fourier space
    %t=+1 +K; -1 -K
    w=parameters.w;
    t=parameters.valley;
    if t==1
        omega=parameters.omega;
    else
        omega=conj(parameters.omega);
    end
    if parameters.shift==2
        h1=h1+2/3*t; % the seemingly-wrong sign is because we need to revert it backwards to have a same convention as before
        h2=h2+1/3*t;
    end
    re=w*(iseq(h1,0).*iseq(h2,0)...   %gamma point
        +conj(omega)*iseq(h1,1*t).*iseq(h2,0)...    %g2
        +conj(omega^2)*iseq(h1,1*t).*iseq(h2,1*t));   %g3
    end