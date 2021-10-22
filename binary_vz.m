function binary_vz(nu,Nmax,w,Nk,ep,Vb,d,vz_t_list)
    % vz_t_list=[min,max] should have different chern number
    chern_p_0_list=[];
    chern_m_0_list=[];
    gap_final_0_list={};
    tot_final_0_list={};
    s0_0_list={};
    sx_0_list={};
    sy_0_list={};
    sz_0_list={};
    rmap_x_0_list={};
    rmap_y_0_list={};

    chern_p_1_list=[];
    chern_m_1_list=[];
    gap_final_1_list={};
    tot_final_1_list={};
    s0_1_list={};
    sx_1_list={};
    sy_1_list={};
    sz_1_list={};
    rmap_x_1_list={};
    rmap_y_1_list={};

    epoch_list=vz_t_list*0;
    ave1=0;
    ave2=0;
    epoch=0;
    
    params=mainTMD('Nmax',Nmax,'V_t',0,'psi_t',240,'V_b',Vb,'psi_b',-14,'vz_t',vz_t_list(1),'vz_b',0,'w',w,'nu',nu,'n',Nk,'epsilon',ep,'d',d);
    
    

    for vz_t_index=1:length(vz_t_list)
        vz_t=vz_t_list(vz_t_index);
        fprintf("vz_t=%.1f\n",vz_t);
        ave1=0;
        ave2=0;
        epoch=0;
        [ave1,ave2,gap_list,tot_list,epoch,chern_p,chern_m,s0,sx,sy,sz,rmap_x,rmap_y,chern_cond]=iter(vz_t,ave1,ave2,epoch,params);
        if abs(chern_p)>0.5
            assert(isempty(chern_p_1_list),'both are on the topological side')
            chern_p_1_list(1)=chern_p;
            chern_m_1_list(1)=chern_m;
            gap_final_1_list{1}=gap_list;
            tot_final_1_list{1}=tot_list;
            % innergap_1_list(1)=innergap;
            % epoch_1_list(1)=epoch;
            s0_1_list{1}=s0;
            sx_1_list{1}=sx;
            sy_1_list{1}=sy;
            sz_1_list{1}=sz;
            rmap_x_1_list{1}=rmap_x;
            rmap_y_1_list{1}=rmap_y;
            ave1_1=ave1;
            ave2_1=ave2;
            epoch_1=epoch;
        else
            assert(isempty(chern_p_0_list),'both are on the trivial side')
            chern_p_0_list(1)=chern_p;
            chern_m_0_list(1)=chern_m;
            gap_final_0_list{1}=gap_list;
            tot_final_0_list{1}=tot_list;
            % innergap_0_list(1)=innergap;
            % epoch_0_list(1)=epoch;
            s0_0_list{1}=s0;
            sx_0_list{1}=sx;
            sy_0_list{1}=sy;
            sz_0_list{1}=sz;
            rmap_x_0_list{1}=rmap_x;
            rmap_y_0_list{1}=rmap_y;
            ave1_0=ave1;
            ave2_0=ave2;
            epoch_0=epoch;
        end
    end
    left=vz_t_list(1); %assume it's trivial
    right=vz_t_list(2); %assme it's topological
    
    while gap_list(end)>1e-3 & abs(right-left)<1e-6
        middle=(left+right)/2;
        vz_t=middle;
        fprintf("vz_t=%e\n",vz_t);
        % ave1=0;
        % ave2=0;
        [ave1,ave2,gap_list,tot_list,epoch,chern_p,chern_m,s0,sx,sy,sz,rmap_x,rmap_y,chern_cond]=iter(vz_t,ave1_1,ave2_1,epoch,params);
        if abs(chern_p)>.5
            chern_p_1_list(end+1)=chern_p;
            chern_m_1_list(end+1)=chern_m;
            gap_final_1_list{end+1}=gap_list;
            tot_final_1_list{end+1}=tot_list;
            % innergap_1_list(end+1)=innergap;
            % epoch_1_list(end+1)=epoch;
            s0_1_list{end+1}=s0;
            sx_1_list{end+1}=sx;
            sy_1_list{end+1}=sy;
            sz_1_list{end+1}=sz;
            rmap_x_1_list{end+1}=rmap_x;
            rmap_y_1_list{end+1}=rmap_y;
            % epoch_1=epoch;
            ave1_1=ave1;
            ave2_1=ave2;
            right=middle;
        else
            chern_p_0_list(end+1)=chern_p;
            chern_m_0_list(end+1)=chern_m;
            gap_final_0_list{end+1}=gap_list;
            tot_final_0_list{end+1}=tot_list;
            % innergap_0_list(end+1)=innergap;
            % epoch_0_list(end+1)=epoch;
            s0_0_list{end+1}=s0;
            sx_0_list{end+1}=sx;
            sy_0_list{end+1}=sy;
            sz_0_list{end+1}=sz;
            rmap_x_0_list{end+1}=rmap_x;
            rmap_y_0_list{end+1}=rmap_y;
            % epoch_0=epoch;
            ave1_0=ave1;
            ave2_0=ave2;
            left=middle;
        end
    end



    save(sprintf('phase_nu%d,%d_ep%.1f_w%.1f_Vb%.1f_d%.1f_Nk%d.mat',nu(1),nu(2),ep,w,Vb,d,Nk),'chern_p_0_list','chern_m_0_list','gap_final_0_list','tot_final_0_list','vz_t_0_list','s0_0_list','sx_0_list','sy_0_list','sz_0_list','rmap_x_0_list','rmap_y_0_list','chern_p_1_list','chern_m_1_list','gap_final_1_list','tot_final_1_list','vz_t_1_list','s0_1_list','sx_1_list','sy_1_list','sz_1_list','rmap_x_1_list','rmap_y_1_list');
end



function [ave1,ave2,gap_list,tot_list,epoch,chern_p,chern_m,s0,sx,sy,sz,rmap_x,rmap_y,chern_cond]=iter(vz_t,ave1,ave2,epoch,params)
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
    for i=1:1000
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
    [~,~,fig_band,fig_spin,chern_p,chern_m,s0,sx,sy,sz,rmap_x,rmap_y,chern_cond]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,strcat('fsgx',params.chern),params);
    savefig(fig_band,strcat(fn,'_band.fig'));
    savefig(fig_spin,strcat(fn,'_spin.fig'));
end
    