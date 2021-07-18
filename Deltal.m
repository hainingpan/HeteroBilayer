function re=Deltal(h1,h2,l,parameters)
    %l=1 for bottom layer
    %l=-1 for top layer
    %t=1 for +K
    %t=-1 for -K
    if l==1
        V=parameters.V_b;
        psi=parameters.psi_b;
    else
        V=parameters.V_t;
        psi=parameters.psi_t;
    end
    t=parameters.valley;
    re=V*((h1==0).*(h2==-1*t)*exp(1i*t*psi)+...    %G1
        (h1==0).*(h2==1*t)*exp(-1i*t*psi)+...    %G4
        (h1==-1*t).*(h2==0)*exp(1i*t*psi)+...  `   %G5
        (h1==1*t).*(h2==0)*exp(-1i*t*psi)+...    %G2
        (h1==1*t).*(h2==1*t)*exp(1i*t*psi)+...    %G3
        (h1==-1*t).*(h2==-1*t)*exp(-1i*t*psi));       %G6
    end