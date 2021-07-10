function [ldos,enlist]=LDOS_TMD_rx(parameters)
    %AA:r=(0,0)
    %AB:r=(1/sqrt(3) aM,0) aM=a/theta
    % rlistx=linspace(0,sqrt(3),100)*parameters.a/parameters.theta;
    % rlisty=0*rlistx;
    r=[0,0];
    rAA=(parameters.aM1+parameters.aM2);
    rlistx=linspace(r(1),rAA(1),40);
    rlisty=linspace(r(2),rAA(2),40);
    eta=1e-3;
    
    n=20;
    xrange=-n:n;
    yrange=-n:n;
    bM1=parameters.bM1;
    bM2=parameters.bM2;
    % kp=parameters.kp;
    % kn=parameters.kn;
    a1=-bM1/(2*n);
    a2=(bM1+bM2)/(2*n);
    enmap=zeros(2*n+1,2*n+1,2*(2*parameters.Nmax+1)^2);
    
    Nrx=length(rlistx);
    Nx=length(xrange);
    Ny=length(yrange);
    V=abs(cross([parameters.aM1,0],[parameters.aM2,0]));
    V=V(3);
    
    psir=zeros(Nx,Ny,length(rlistx),2*(2*parameters.Nmax+1)^2);
    for xindex=1:Nx
        kx=xrange(xindex);
        for yindex=1:Ny
            ky=yrange(yindex);
            k=kx*a1+ky*a2;
            [val,vec]=energyTMD(k(1),k(2),parameters);
            enmap(xindex,yindex,:)=val;
    %         for i=1:length(rlistx)
    %             psir{i}(xindex,yindex,:)=u2(vec,[rlistx(i),0],parameters);
    %         end
            [psib,psit]=u(vec,rlistx,rlisty,parameters);
            psib=squeeze(psib)/(sqrt(V)/5.076e-3);
            psit=squeeze(psit)/(sqrt(V)/5.076e-3);
            psir(xindex,yindex,:,:)=abs(psib).^2+abs(psit).^2;
        end
    end
    
    % enlist=linspace(min(enmap(:,:,4),[],'all'),1.2*max(enmap(:,:,1),[],'all'),100);
    enlist=linspace(-60e-3,0e-3,100);
    Nen=length(enlist);
    ldos=zeros(length(rlistx),length(enlist));
    
    parfor i=1:Nrx
        for j=1:Nen
            deltaf=1/pi*eta./((enlist(j)-enmap).^2+eta^2);
            ldosprod=squeeze(psir(:,:,i,:)).*deltaf;
            ldos(i,j)=sum(ldosprod(:)); %in the unit of eV^-1*nm^-2
        end
    end
    
    ldos=ldos/(Nx*Ny);
    
    end
    
    