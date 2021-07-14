function re=DeltaTT(h1,h2,parameters)
    %Delta_T^dagger in Fourier space
    w=parameters.w;
    omega=parameters.omega;
    re=w*((h1==0).*(h2==0)...   %gamma point
        +conj(omega)*(h1==-1).*(h2==0)...    %g2
        +conj(omega^2)*(h1==-1).*(h2==-1));   %g3
    end