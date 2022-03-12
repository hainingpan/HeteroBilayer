% function iter_SDW(nu,Nmax)
Nmax=2;
nu=[2,1];
w=12;
vz=-20;
Vb=7;
ep=2.5e2;

params=mainTMD('SDW',10e-3,'d',5,'Nmax',Nmax,'V_t',0,'psi_t',240,'V_b',Vb,'psi_b',-14,'vz_t',vz,'vz_b',0,'w',w,'nu',nu,'n',15,'epsilon',ep,'shift',2);
[energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(0,0,0,params);
[ave1,ave2,occ]=average(energyall,wfall,0,params); 

[gap,tot,fig1,fig2]=plotline_2(energyall,0,0,V1_ave_delta,V2_ave_delta,ave1,ave2,0,strcat('fs',params.chern),params);
gap_list=[gap];
tot_list=[tot];
fprintf("%d: Gap=%e (meV), E=%e (meV)\n",0,gap*1e3,tot*1e3);
fig_tot=figure;
for i=1:100
    [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,i,params);    
    [ave1_n,ave2_n,occ]=average(energyall,wfall,i,params);
    [gap,tot]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,i,strcat('fs',params.chern),params);
    gap_list(end+1)=gap;
    tot_list(end+1)=tot;
    pct_change=((tot_list(end)-tot_list(end-1))/tot_list(end-1));
    fprintf("%d: Gap=%.10f (meV), E=%.10f (meV), change=%e\n",i,gap*1e3,tot*1e3,pct_change);
    figure(fig_tot);
    plot(tot_list(2:end)*1e3);
    ave1=ave1_n;
    ave2=ave2_n;
    if abs(pct_change)<1e-7
        break;
    end
end

