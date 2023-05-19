function [ldos,enlist,rlist]=LDOS_TMD_rx_wf(energyall,wfall,parameters)
    %AA:r=(0,0)
    %AB:r=(1/sqrt(3) aM,0) aM=a/theta
    % rlistx=linspace(0,sqrt(3),100)*parameters.a/parameters.theta;
    % rlisty=0*rlistx;
    r=[0,0];
    rAA=(parameters.aM1+parameters.aM2);
    rlistx=linspace(r(1),rAA(1),40);
    rlisty=linspace(r(2),rAA(2),40);
    rlist=sqrt(rlistx.^2+rlisty.^2);
    eta=1e-3;
    Nk=size(parameters.k,1);
    
    n=20;
    xrange=-n:n-1;
    yrange=-n:n-1;
    bM1=parameters.bM1;
    bM2=parameters.bM2;
%     a1=bM1/(2*n);
%     a2=bM2/(2*n);
    a1=-bM1/(2*n);
    a2=(bM1+bM2)/(2*n);
    % enmap=zeros(Nk,4*size(parameters.h1index,1));
    enmap=energyall;
    
    Nrx=length(rlistx);
    Nx=length(xrange);
    Ny=length(yrange);
    V=abs(cross([parameters.aM1,0],[parameters.aM2,0]));
    V=V(3);
    
    psir=zeros(Nk,length(rlistx),4*size(parameters.h1index,1));
    for k_index=1:Nk
        vec=squeeze(wfall(k_index,:,:)).';
        [psi_b_p,psi_t_p,psi_b_m,psi_t_m]=u_wf(vec,rlistx,rlisty,parameters);
        psi_b_p=squeeze(psi_b_p)/(sqrt(V)/5.076e-3);
        psi_t_p=squeeze(psi_t_p)/(sqrt(V)/5.076e-3);
        psi_b_m=squeeze(psi_b_m)/(sqrt(V)/5.076e-3);
        psi_t_m=squeeze(psi_t_m)/(sqrt(V)/5.076e-3);
        psir(k_index,:,:)=abs(psi_b_p).^2+abs(psi_t_p).^2+abs(psi_b_m).^2+abs(psi_t_m).^2;
        % psir(k_index,:,:)=abs(psi_b_p).^2+abs(psi_t_p).^2;
    end
    
    % enlist=linspace(min(enmap(:,5),[],'all')-10e-3,10e-3+max(enmap(:,1),[],'all'),200);
    enlist=linspace(min(enmap(:,1),[],'all')-5e-3,10e-3+max(enmap(:,5),[],'all'),200);
%     enlist=linspace(40e-3,50e-3,100);
    Nen=length(enlist);
    ldos=zeros(length(rlistx),length(enlist));
    
    for i=1:Nrx
        for j=1:Nen
            deltaf=1/pi*eta./((enlist(j)-enmap).^2+eta^2);
            ldosprod=squeeze(psir(:,i,:)).*deltaf;
            ldos(i,j)=sum(ldosprod(:)); %in the unit of eV^-1*nm^-2
        end
    end
    
    ldos=ldos/(Nx*Ny);
    
    end
    
    