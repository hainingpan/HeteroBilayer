function [ldos,enlist,rlist]=LDOS_uc(energyall,wfall,params)
    r=[0,0];
    rAA=(params.aM1+params.aM2);
    rlistx=linspace(r(1),rAA(1),40);
    rlisty=linspace(r(2),rAA(2),40);
    rlist=sqrt(rlistx.^2+rlisty.^2);
    eta=2e-3;
    Nk=size(params.k,1);
    Nb=size(params.b,1);
    Nrx=size(rlistx,2);
    NL=size(wfall,2);
    V=params.area;

    psir=zeros(Nk,Nrx,NL);
    for k_index=1:Nk
        expo=exp(1i*params.b*[rlistx(:),rlisty(:)]');
        vec_b_p=squeeze(wfall(k_index,:,1:Nb));
        vec_t_p=squeeze(wfall(k_index,:,Nb+1:2*Nb));
        vec_b_m=squeeze(wfall(k_index,:,2*Nb+1:3*Nb));
        vec_t_m=squeeze(wfall(k_index,:,3*Nb+1:4*Nb));
        u_b_p=vec_b_p*expo;
        u_t_p=vec_t_p*expo;
        u_b_m=vec_b_m*expo;
        u_t_m=vec_t_m*expo;
        psir(k_index,:,:)=abs(u_b_p').^2+abs(u_t_p').^2+abs(u_b_m').^2+abs(u_t_m').^2;
    end

    enlist=linspace(min(energyall(:,1)),max(energyall(:,2)),100);
    Nen=length(enlist);
    ldos=zeros(length(rlistx),length(enlist));

    for i=1:Nrx
        for j=1:Nen
            deltaf=1/pi*eta./((enlist(j)-energyall).^2+eta.^2);
            ldosprod=squeeze(psir(:,i,:)).*deltaf;
            ldos(i,j)=sum(ldosprod(:));
        end
    end

    ldos=ldos/(Nk);
end
