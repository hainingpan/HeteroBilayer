function re=DeltaT(h1,h2,parameters)
    %Delta_T in Fourier space
    w=parameters.w;
    omega=parameters.omega;
    re=w*((h1==0).*(h2==0)+...  %gamma point
        omega*(h1==1).*(h2==0)...  %g2
        +omega^2*(h1==1).*(h2==1));   %g3
    end