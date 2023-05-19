function sweep_ep(nu,Nmax,w,Nk,Vz_t,Vb,d)
    ep_list=linspace(150,15,101);
    % ep_list=1./linspace(0,1/50,20);
    chern_p_list=ep_list;
    chern_m_list=ep_list;
    gap_final_list={};
    tot_final_list={};
    s0_list={};
    sx_list={};
    sy_list={};
    sz_list={};
    rmap_x_list={};
    rmap_y_list={};
    % chern_cond_list=ep_list;
   pol_list=ep_list;

    epoch_list=ep_list;
    ave1=0;
    ave2=0;
    epoch=0;
    
    params=mainTMD('Nmax',Nmax,'V_t',0,'psi_t',240,'V_b',Vb,'psi_b',-14,'vz_b',0,'w',w,'nu',nu,'n',Nk,'Vz_t',Vz_t,'d',d);
    for ep_index=1:length(ep_list)
        epsilon=ep_list(ep_index);
        fprintf("epsilon=%.1f\n",epsilon);
%         ave1=0;
%         ave2=0;
%         epoch=0;
        [ave1,ave2,gap_list,tot_list,epoch,chern_p,chern_m,s0,sx,sy,sz,rmap_x,rmap_y,chern_cond,pol]=iter(epsilon,ave1,ave2,epoch,params);

        chern_p_list(ep_index)=chern_p;
        chern_m_list(ep_index)=chern_m;
        gap_final_list{ep_index}=gap_list;
        tot_final_list{ep_index}=tot_list;
        % innergap_list(ep_index)=innergap;
        epoch_list(ep_index)=epoch;
        s0_list{ep_index}=s0;
        sx_list{ep_index}=sx;
        sy_list{ep_index}=sy;
        sz_list{ep_index}=sz;
        rmap_x_list{ep_index}=rmap_x;
        rmap_y_list{ep_index}=rmap_y;
        % chern_cond_list(ep_index)=chern_cond;
        pol_list(ep_index)=pol;
    end

    save(sprintf('phase_nu%d,%d_Vz%.1f_w%.1f_Vb%.1f_d%.1f_Nk%d.mat',nu(1),nu(2),Vz_t,w,Vb,d,Nk),'chern_p_list','chern_m_list','gap_final_list','tot_final_list','epoch_list','ep_list','s0_list','sx_list','sy_list','sz_list','rmap_x_list','rmap_y_list','pol_list');
end



function [ave1,ave2,gap_list,tot_list,epoch,chern_p,chern_m,s0,sx,sy,sz,rmap_x,rmap_y,chern_cond,pol]=iter(epsilon,ave1,ave2,epoch,params)
    params.epsilon=epsilon;
    [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,epoch,params);
    [ave1_n,ave2_n,occ]=average(energyall,wfall,epoch,params); 
    
    re=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,0,'',params);
    gap=re.gap;
    tot=re.tot;
    
    gap_list=[gap];
    tot_list=[tot];
    fprintf("%d: Gap=%.10f (meV), E=%.10f (meV)\n",0,gap*1e3,tot*1e3);
    ave1=ave1_n;
    ave2=ave2_n;    
    epoch=epoch+1;
    for i=1:1000
        [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,epoch,params);
        [ave1_n,ave2_n,occ]=average(energyall,wfall,epoch,params);       
        re=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,'',params);
        gap=re.gap;
        tot=re.tot;
        gap_list(end+1)=gap;
        tot_list(end+1)=tot;
        ave1=ave1_n;
        ave2=ave2_n;
        drawnow;
        pct_change=((tot_list(end)-tot_list(end-1))/tot_list(end-1));
        fprintf("%d: Gap=%.10f (meV), E=%.10f (meV), change=%e\n",i,gap*1e3,tot*1e3,pct_change);
        if abs(pct_change)<1e-9
            break;
        end
        epoch=epoch+1;
    end
    fn=sprintf('nu_%d,%d_Nmax%d_w%.1f_Vb%.1f_Nk_%d_Vzt_%.1f_d%.1f_ep%.2f',params.nu(1),params.nu(2),params.Nmax,1000*params.w,1000*params.V_b,params.n,params.Vz_t*1000,params.d/(1e-9*5.076e6),params.epsilon);
    re=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,strcat('fgsp',params.chern),params);
    gap=re.gap;
    tot=re.tot;
    fig_band=re.fig_band;
    fig_spin=re.fig_spin;
    chern_p=re.chern_p;
    chern_m=re.chern_m;
    s0=re.s0_list;
    sx=re.sx_list;
    sy=re.sy_list;
    sz=re.sz_list;
    rmap_x=re.rmap_x;
    rmap_y=re.rmap_y;
    chern_cond=re.chern_cond;
    pol=re.pol;
    gap_list(end+1)=gap;
    tot_list(end+1)=tot;
    savefig(fig_band,strcat(fn,'_band.fig'));
    savefig(fig_spin,strcat(fn,'_spin.fig'));
end
    