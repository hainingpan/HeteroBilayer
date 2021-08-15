function sweep_vz(nu,Nmax,w,Nk,ep,vz_t_list)
    % vz_t_list=-50:50;
    chern_p_list=vz_t_list;
    chern_m_list=vz_t_list;
    gap_final_list={};
    tot_final_list={};
    epoch_list=vz_t_list;
    ave1=0;
    ave2=0;
    epoch=0;

    for vz_t_index=1:length(vz_t_list)
        vz_t=vz_t_list(vz_t_index);
        fprintf("vz_t=%.1f\n",vz_t);
        [ave1,ave2,gap_list,tot_list,epoch,chern_p,chern_m]=iter(nu,Nmax,w,Nk,vz_t,ep,ave1,ave2,epoch);

        chern_p_list(vz_t_index)=chern_p;
        chern_m_list(vz_t_index)=chern_m;
        gap_final_list{vz_t_index}=gap_list;
        tot_final_list{vz_t_index}=tot_list;
        % innergap_list(vz_t_index)=innergap;
        epoch_list(vz_t_index)=epoch;
    end

    save(sprintf('phase_nu%d,%d_ep%d_w%.1f_Nk_%d.mat',params.nu(1),params.nu(2),epsilon,w,Nk),'chern_p_list','chern_m_list','gap_final_list','tot_final_list','epoch_list','vz_t_list')




function [ave1,ave2,gap_list,tot_list,epoch,chern_p,chern_m]=iter(nu,Nmax,w,Nk,vz_t,ep,ave1,ave2,epoch)
    params=mainTMD('SDW',20e-3,'Nmax',Nmax,'V_t',0,'psi_t',240,'V_b',15,'psi_b',-14,'vz_t',vz_t,'vz_b',0,'w',w,'nu',nu,'n',Nk,'epsilon',ep,'shift',1);
    [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,epoch,params);
    [ave1_n,ave2_n,occ]=average(energyall,wfall,epoch,params); 
    
    [gap,tot]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,0,'',params);
    
    gap_list=[gap];
    tot_list=[tot];
    fprintf("%d: Gap=%.10f (meV), E=%.10f (meV)\n",0,gap*1e3,tot*1e3);
    ave1=ave1_n;
    ave2=ave2_n;    
    epoch=epoch+1;
    for i=1:500
        [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,epoch,params);
        [ave1_n,ave2_n,occ]=average(energyall,wfall,epoch,params);       
        [gap,tot]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,'',params);
        gap_list(end+1)=gap;
        tot_list(end+1)=tot;
        ave1=ave1_n;
        ave2=ave2_n;
        drawnow;
        pct_change=((tot_list(end)-tot_list(end-1))/tot_list(end-1));
        fprintf("%d: Gap=%.10f (meV), E=%.10f (meV), change=%e\n",i,gap*1e3,tot*1e3,pct_change);
        if abs(pct_change)<1e-7
            break;
        end
        epoch=epoch+1;
    end
    fn=sprintf('nu_%d,%d_Nmax%d_w%.1f_Nk_%d_Vzt_%.1f_ep%.1f',params.nu(1),params.nu(2),params.Nmax,1000*params.w,Nk,vz_t,ep);
    [~,~,fig_band,fig_spin,chern_p,chern_m]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,strcat('fsg',params.chern),params);
    savefig(fig_band,strcat(fn,'_band.fig'));
    savefig(fig_spin,strcat(fn,'_spin.fig'));
    % save(strcat(fn,'.mat'),'gap_list','tot_list','i','chern_p','chern_m');
    
    