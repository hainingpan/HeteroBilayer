% function iter_SDW(nu,Nmax)
Nmax=3;
nu=[6,6];
params=mainTMD('Nmax',Nmax,'V_t',0,'psi_t',240,'V_b',15,'psi_b',-14,'vz_t',-4,'vz_b',0,'w',20,'nu',nu,'n',27,'epsilon',25,'shift',1);
[energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(0,0,0,params);
[ave1,ave2,occ]=average(energyall,wfall,0,params); 

[gap,tot,fig1,fig2]=plotline_2(energyall,0,0,V1_ave_delta,V2_ave_delta,ave1,ave2,0,strcat('fs',params.chern),params);
fig_tot=figure;
gap_list=[gap];
tot_list=[tot];
fprintf("%d: Gap=%e (meV), E=%e (meV)\n",0,gap*1e3,tot*1e3);
for i=1:100
    [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,i,params);    
    [ave1_n,ave2_n,occ]=average(energyall,wfall,i,params);
    [gap,tot]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,i,strcat('fs',params.chern),params);
    gap_list(end+1)=gap;
    tot_list(end+1)=tot;
    fprintf("%d: Gap=%e (meV), E=%e (meV)\n",i,gap*1e3,tot*1e3);
    figure(fig_tot);
    plot(tot_list(2:end)*1e3);
    ave1=ave1_n;
    ave2=ave2_n;
    if abs(tot_list(end)-tot_list(end-1))<1e-10
        break;
    end
end

