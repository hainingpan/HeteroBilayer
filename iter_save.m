function iter_save(nu,Nmax)
    params=mainTMD('Nmax',Nmax,'V_t',0,'psi_t',240,'V_b',15,'psi_b',-14,'vz_t',-40,'vz_b',0,'w',0,'nu',nu,'n',15,'epsilon',25,'shift',1);
    [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(0,0,0,params);
    [ave1,ave2,occ]=average(energyall,wfall,0,params); 
    [X_index,Y_index]=meshgrid(linspace(0,1,22));
    [sitesX_index,sitesY_index]=meshgrid(-1:3);
    rmap=[X_index(:),Y_index(:)]*[params.am1;params.am2];
    rsite=[sitesX_index(:),sitesY_index(:)]*[params.aM1;params.aM2];
    rmap_x=rmap(:,1);
    rmap_y=rmap(:,2);
    
    [gap,tot]=plotline_2(energyall,0,0,V1_ave_delta,V2_ave_delta,ave1,ave2,0,'',params);
    
    gap_list=[gap];
    tot_list=[tot];
    drawnow;
    for i=1:100
        disp(i)
        [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,i,params);    
        [ave1_n,ave2_n,occ]=average(energyall,wfall,i,params);
        
        [gap,tot]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,i,'',params);
        gap_list(end+1)=gap;
        tot_list(end+1)=tot;
        ave1=ave1_n;
        ave2=ave2_n;
        drawnow;
        if abs(tot_list(end)-tot_list(end-1))<1e-10
            break;
        end
    end
    
    [gap,tot,fig_band]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,i,'f',params);
    savefig(fig_band,sprintf('nu_%d,%d_Nmax%d_band.fig',params.nu(1),params.nu(2),params.Nmax))
    fig1=figure;
    [s0,sx,sy,sz]=S_r(ave2_n,rmap_x,rmap_y,1,params);
    quiver(rmap_x(:),rmap_y(:),(sx(:)),(sy(:)));
    daspect([1,1,1]);
    hold on;
    scatter(rsite(:,1),rsite(:,2));
    xlim([min(rmap_x)*1.1,max(rmap_x)*1.1]);
    ylim([min(rmap_y)*1.1,max(rmap_y)*1.1]);
    savefig(fig1,sprintf('nu_%d,%d_Nmax%d_spin.fig',params.nu(1),params.nu(2),params.Nmax))

    save(sprintf('nu_%d,%d_Nmax%d.mat',params.nu(1),params.nu(2),params.Nmax),'gap_list','tot_list','i');
    
    