function [ubgridp,utgridp]=u(vec,rx,ry,parameters)
    %rx,ry: scalar|array real space position
    [~,Nstate]=size(vec);
    [Nx,Ny]=size(rx);
    Nmax=parameters.Nmax;
    Nrange=-Nmax:Nmax;
    % [h1index,h2index]=meshgrid(Nrange,Nrange);
    h1index=parameters.h1index;
    h2index=parameters.h2index;

    
    % V=abs(cross([parameters.aM1,0],[parameters.aM2,0]));
    % V=V(3);
    
    expo=exp(1i*(h1index(:)*parameters.bM1+h2index(:)*parameters.bM2)*[rx(:),ry(:)]');
    vecb=vec(1:size(parameters.h1index,1),:);
    vect=vec(size(parameters.h1index,1)+1:2*size(parameters.h1index,1),:);
    ub=vecb'*expo;
    ut=vect'*expo;
    
    ubgrid=reshape(ub,Nstate,Nx,Ny);
    utgrid=reshape(ut,Nstate,Nx,Ny);    %Normalization to V_UC
    
    ubgridp=permute(ubgrid,[2,3,1]);
    utgridp=permute(utgrid,[2,3,1]);
    
    
    end