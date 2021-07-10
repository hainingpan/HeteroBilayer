function [ubgridp,utgridp]=u(vec,rx,ry,parameters)
    %rx,ry: scalar|array real space position
    [~,Nstate]=size(vec);
    [Nx,Ny]=size(rx);
    Nmax=parameters.Nmax;
    Nrange=-Nmax:Nmax;
    [h1index,h2index]=meshgrid(Nrange,Nrange);
    
    % V=abs(cross([parameters.aM1,0],[parameters.aM2,0]));
    % V=V(3);
    
    expo=exp(1i*(h1index(:)*parameters.bM1+h2index(:)*parameters.bM2)*[rx(:),ry(:)]');
    vecb=vec(1:(2*Nmax+1)^2,:);
    vect=vec((2*Nmax+1)^2+1:2*(2*Nmax+1)^2,:);
    ub=vecb'*expo;
    ut=vect'*expo;
    
    ubgrid=reshape(ub,Nstate,Nx,Ny);
    utgrid=reshape(ut,Nstate,Nx,Ny);    %Normalization to V_UC
    
    ubgridp=permute(ubgrid,[2,3,1]);
    utgridp=permute(utgrid,[2,3,1]);
    
    
    end