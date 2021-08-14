function iter_save(nu,Nmax,w,Nk,vz_t,SDW)
    params=mainTMD('SDW',SDW,'Nmax',Nmax,'V_t',0,'psi_t',240,'V_b',15,'psi_b',-14,'vz_t',vz_t,'vz_b',0,'w',w,'nu',nu,'n',Nk,'epsilon',25,'shift',1);
    [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(0,0,0,params);
    [ave1,ave2,occ]=average(energyall,wfall,0,params); 
    
    [gap,tot]=plotline_2(energyall,0,0,V1_ave_delta,V2_ave_delta,ave1,ave2,0,'',params);
    
    gap_list=[gap];
    tot_list=[tot];
    fprintf("%d: Gap=%e (meV), E=%e (meV)",0,gap,tot);
    for i=1:500
        [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,i,params);    
        [ave1_n,ave2_n,occ]=average(energyall,wfall,i,params);
        
        [gap,tot]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,i,'',params);
        gap_list(end+1)=gap;
        tot_list(end+1)=tot;
        fprintf("%d: Gap=%e (meV), E=%e (meV)",i,gap,tot);
        ave1=ave1_n;
        ave2=ave2_n;
        drawnow;
        if abs(tot_list(end)-tot_list(end-1))<1e-10
            break;
        end
    end
    fn=sprintf('nu_%d,%d_Nmax%d_w%.1f_Nk_%d_Vzt_%.1f',params.nu(1),params.nu(2),params.Nmax,1000*params.w,Nk,vz_t);
    [gap,tot,fig_band,fig_spin]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,i,strcat('fs',params.chern),params);
    savefig(fig_band,strcat(fn,'_band.fig'))
    savefig(fig_spin,strcat(fn,'_spin.fig'))

    save(strcat(fn,'.mat'),'gap_list','tot_list','i');
    
    