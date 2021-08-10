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
    re=w*(iseq(h1,0).*iseq(h2,0)...   %gamma point
        +conj(omega)*iseq(h1,1*t).*iseq(h2,0)...    %g2
        +conj(omega^2)*iseq(h1,1*t).*iseq(h2,1*t));   %g3
    end