% params=mainTMD('Nmax',3,'d',100e-9*5.076e6,'a_t',3.42e-10*5.076e6,'a_b',3.575e-10*5.076e6,'V_t',20.1,'psi_t',240,'V_b',20.1,'psi_b',-14,'vz_t',-100,'vz_b',0,'w',0,'nu',[1,1],'n',15,'epsilon',20);
% params=mainTMD('Nmax',3,'V_t',0,'psi_t',240,'V_b',15,'psi_b',-14,'vz_t',-4,'vz_b',0,'w',20,'nu',[2,1],'n',15,'epsilon',20,'tsymm',1);
params=mainTMD('Nmax',3,'V_t',0,'psi_t',240,'V_b',15,'psi_b',-14,'vz_t',-40,'vz_b',0,'w',0,'nu',[2,2],'n',15,'epsilon',25,'shift',1);
[energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(0,0,0,params);

[ave1,ave2,occ]=average(energyall,wfall,0,params); 
% [ldos,enlist,rlist]=LDOS_uc(energyall,wfall,params);
% figure;surf(rlist/params.aM,1000*enlist,ldos','edgecolor','none');view(2);
fig1=figure;
% gap=plotline(energyall,wfall,occ,valley_index,0,0,params);
[gap,tot]=plotline_2(energyall,0,0,V1_ave_delta,V2_ave_delta,ave1,ave2,0,params);

gap_list=[gap];
tot_list=[tot];


for i=1:100
    disp(i)
    [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,i,params);    
    [ave1_n,ave2_n,occ]=average(energyall,wfall,i,params);
    % gap=plotline(energyall,wfall,occ,valley_index,ave1,ave2,params);
    [gap,tot]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,i,params);
    gap_list(end+1)=gap;
    tot_list(end+1)=tot;
    figure(fig1);
    plot(tot_list(2:end)*1e3);
    ave1=ave1_n;
    ave2=ave2_n;
    if abs(tot_list(end)-tot_list(end-1))<1e-8
        break;
    end
end

save('final_2,2.mat','gap_list','tot_list','energyall','ave1','ave2','V1_ave_delta','V2_ave_delta','ave1_n','ave2_n','i','params');


