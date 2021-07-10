function re=DeltaT(h1,h2,parameters)
    %Delta_T in Fourier space
    w=parameters.w;
    re=w*((h1==0).*(h2==0)+exp(1i*2*pi/3)*(h1==-1).*(h2==0)+exp(1i*4*pi/3)*(h1==-1).*(h2==-1));
    end