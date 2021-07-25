% params=mainTMD('Nmax',3,'d',100e-9*5.076e6,'a_t',3.42e-10*5.076e6,'a_b',3.575e-10*5.076e6,'V_t',20.1,'psi_t',240,'V_b',20.1,'psi_b',-14,'vz_t',-100,'vz_b',0,'w',0,'nu',[1,1],'n',15,'epsilon',20);
params=mainTMD('Nmax',3,'V_t',16.1,'psi_t',240,'V_b',16.1,'psi_b',-14,'vz_t',-10,'vz_b',0,'w',8.3,'nu',[2,1],'n',15,'epsilon',20);
[energyall,wfall,valley_index]=energyMF(0,0,params);
% [ave1,ave2,occ,~]=average(energyall,wfall,1,0,params); %for QAHE

[ave1,ave2,occ]=average(energyall,wfall,0,1,params); %for QSHE
% [ldos,enlist,rlist]=LDOS_uc(energyall,wfall,params);
% figure;surf(rlist/params.aM,1000*enlist,ldos','edgecolor','none');view(2);
% fig1=figure;
gap=plotline(energyall,wfall,occ,valley_index,0,0,params);

% % gap_list=[gap];


for i=1:100
    [energyall,wfall,valley_index]=energyMF(ave1,ave2,params);    
    [ave1_n,ave2_n,occ]=average(energyall,wfall,0,1,params);
    gap=plotline(energyall,wfall,occ,valley_index,ave1,ave2,params);
    gap_list(end+1)=gap;
%     figure(fig1);
%     plot(gap_list*1e3);
    ave1=ave1_n;
    ave2=ave2_n;
end

