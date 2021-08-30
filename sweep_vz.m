function sweep_vz(nu,Nmax,w,Nk,ep,Vb,d,vz_t_list)
    % vz_t_list=-50:50;
    chern_p_list=vz_t_list;
    chern_m_list=vz_t_list;
    gap_final_list={};
    tot_final_list={};
    s0_list={};
    sx_list={};
    sy_list={};
    sz_list={};
    rmap_x_list={};
    rmap_y_list={};

    epoch_list=vz_t_list;
    ave1=0;
    ave2=0;
    epoch=0;
    
    params=mainTMD('Nmax',Nmax,'V_t',0,'psi_t',240,'V_b',Vb,'psi_b',-14,'vz_t',vz_t_list(1),'vz_b',0,'w',w,'nu',nu,'n',Nk,'epsilon',ep,'d',d);
    for vz_t_index=1:length(vz_t_list)
        vz_t=vz_t_list(vz_t_index);
        fprintf("vz_t=%.1f\n",vz_t);
%         ave1=0;
%         ave2=0;
%         epoch=0;
        [ave1,ave2,gap_list,tot_list,epoch,chern_p,chern_m,s0,sx,sy,sz,rmap_x,rmap_y]=iter(vz_t,ave1,ave2,epoch,params);

        chern_p_list(vz_t_index)=chern_p;
        chern_m_list(vz_t_index)=chern_m;
        gap_final_list{vz_t_index}=gap_list;
        tot_final_list{vz_t_index}=tot_list;
        % innergap_list(vz_t_index)=innergap;
        epoch_list(vz_t_index)=epoch;
        s0_list{vz_t_index}=s0;
        sx_list{vz_t_index}=sx;
        sy_list{vz_t_index}=sy;
        sz_list{vz_t_index}=sz;
        rmap_x_list{vz_t_index}=rmap_x;
        rmap_y_list{vz_t_index}=rmap_y;
    end

    save(sprintf('phase_nu%d,%d_ep%d_w%.1f_Vb%.1f_d%.1f_Nk%d.mat',nu(1),nu(2),ep,w,Vb,d,Nk),'chern_p_list','chern_m_list','gap_final_list','tot_final_list','epoch_list','vz_t_list','s0_list','sx_list','sy_list','sz_list','rmap_x_list','rmap_y_list')
end



function [ave1,ave2,gap_list,tot_list,epoch,chern_p,chern_m,s0,sx,sy,sz,rmap_x,rmap_y]=iter(vz_t,ave1,ave2,epoch,params)
    params.Vz_t=vz_t*1e-3;
    [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,epoch,params);
    [ave1_n,ave2_n,occ]=average(energyall,wfall,epoch,params); 
    
    [gap,tot]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,0,'',params);
    
    gap_list=[gap];
    tot_list=[tot];
    fprintf("%d: Gap=%.10f (meV), E=%.10f (meV)\n",0,gap*1e3,tot*1e3);
    ave1=ave1_n;
    ave2=ave2_n;    
    epoch=epoch+1;
    for i=1:100
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
    fn=sprintf('nu_%d,%d_Nmax%d_w%.1f_Vb%.1f_Nk_%d_Vzt_%.1f_d%.1f_ep%.1f',params.nu(1),params.nu(2),params.Nmax,1000*params.w,1000*params.V_b,params.n,params.Vz_t*1000,params.d/(1e-9*5.076e6),params.epsilon);
    [~,~,fig_band,fig_spin,chern_p,chern_m,s0,sx,sy,sz,rmap_x,rmap_y]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,strcat('fsg',params.chern),params);
    savefig(fig_band,strcat(fn,'_band.fig'));
    savefig(fig_spin,strcat(fn,'_spin.fig'));
end
    