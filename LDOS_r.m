function [ldos,ldosAB,enlist]=LDOS_r(parameters)
    %:r=(0,0)
    %AB:r=(1/sqrt(3) aM,0) aM=a/theta
    r=[0,0];
    rAA=(parameters.aM1+parameters.aM2);
    rx=linspace(r(1),rAA(1),40);
    ry=linspace(r(2),rAA(2),40);
    eta=1e-3;

    n=10;
    xrange=-n:n;
    yrange=-n:n;
    bM1=parameters.bM1;
    bM2=parameters.bM2;
    a1=bM1/(2*n);
    a2=bM2/(2*n);
    enmap=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);
    psi2=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);
    psiAB2=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);
    
    
    Nx=length(xrange);
    Ny=length(yrange);
    
    for xindex=1:Nx
        kx=xrange(xindex);
        for yindex=1:Ny
            ky=yrange(yindex);
            k=kx*a1+ky*a2;
            [val,vec]=energyTMD(k(1),k(2),parameters);
            enmap(xindex,yindex,:)=val;
            psi2(xindex,yindex,:)=u(vec,rx,ry,parameters);      
        end
    end
    
    enlist=linspace(min(enmap(:,:,4),[],'all'),max(enmap(:,:,1),[],'all'),100);
    ldos=zeros(1,length(enlist));
    for i=1:length(enlist)
    %     fprintf("i_r=%d of %d\n",i,length(enlist));
        deltaf=eta./((enlist(i)-enmap).^2+eta^2);
        ldosprod=psi2.*deltaf;
        ldos(i)=sum(ldosprod(:));
    end 
    ldosAB=zeros(1,length(enlist));
    parfor i=1:length(enlist)
    %     fprintf("i_r=%d of %d\n",i,length(enlist));
        deltaf=eta./((enlist(i)-enmap).^2+eta^2);
        ldosprod=psiAB2.*deltaf;
        ldosAB(i)=sum(ldosprod(:));
    end

    end
    