function re=Deltal(h1,h2,l,parameters)
    %l=1 for bottom layer
    %l=-1 for top layer
    if l==1
        V=parameters.V_b;
        psi=parameters.psi_b;
    else
        V=parameters.V_t;
        psi=parameters.psi_t;
    end
    re=V*((h1==0).*(h2==-1)*exp(1i*psi)+...    %G1
        (h1==0).*(h2==1)*exp(-1i*psi)+...    %G4
        (h1==-1).*(h2==0)*exp(1i*psi)+...  `   %G5
        (h1==1).*(h2==0)*exp(-1i*psi)+...    %G2
        (h1==1).*(h2==1)*exp(1i*psi)+...    %G3
        (h1==-1).*(h2==-1)*exp(-1i*psi));       %G6
    end