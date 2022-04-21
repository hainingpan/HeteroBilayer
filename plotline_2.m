function re=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,output,params)
    re.gap=nan;
    Nai=size(params.ailist,1);  % The expansion of super cell
    % gap
    if ismember('g',output)
        % [energyall_p,energyall_m,~,~]=energyMF_bc(ave1,ave2,'dense',epoch,params);
        [energyall_p,energyall_m,~,~]=energyMF_bc(ave1,ave2,'bc',epoch,params);
        if prod(size(energyall_m))>1
            energyall_p_dense=zeros(size(params.k_dense,1),size(energyall_p,2));
            energyall_m_dense=zeros(size(params.k_dense,1),size(energyall_m,2));
            for level=1:size(energyall_p,2)
                energyall_p_dense(:,level)=griddata(params.k_bc(:,1),params.k_bc(:,2),energyall_p(:,level),params.k_dense(:,1),params.k_dense(:,2));
                energyall_m_dense(:,level)=griddata(params.k_bc(:,1),params.k_bc(:,2),energyall_m(:,level),params.k_dense(:,1),params.k_dense(:,2));
            end
            energyall_sort=sort([energyall_p_dense(:);energyall_m_dense(:)]);
        else
            energyall_p_dense=zeros(size(params.k_dense,1),size(energyall_p,2));
            for level=1:size(energyall_p,2)
                energyall_p_dense(:,level)=griddata(params.k_bc(:,1),params.k_bc(:,2),energyall_p(:,level),params.k_dense(:,1),params.k_dense(:,2));
            end
            energyall_sort=sort(energyall_p_dense(:));
        end
        Nk=size(params.k_dense,1);
    else
        Nk=size(params.k,1);
        energyall_sort=sort(energyall(:));
    end
    mu_v=energyall_sort(Nk*Nai*params.nu(1)/(params.nu(2)));
    mu_c=energyall_sort(1+Nk*Nai*params.nu(1)/(params.nu(2)));
    re.gap=mu_c-mu_v;

    if prod(size(V1_ave_delta))>1 && prod(size(V2_ave_delta))>1
        re.tot=totalenergy(V1_ave_delta,V2_ave_delta,ave1,ave2,epoch,params);
    else
        re.tot=totalenergy(V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,params);
    end

    

    %% Chern number
    if ismember('c',output)
        [re.chern_p,re.chern_m,bcmap_p,bcmap_m]=chern_gs(ave1,ave2,epoch,params);
        chern_str=sprintf('\n(%d)+K:{%.4f} -K:{%.4f}',params.NL,re.chern_p,re.chern_m);
    else
        re.chern_p=nan;
        re.chern_m=nan;
        chern_str='';
    end
    
    %% sigma_xy Hall conductance (problematic, may need to change, suppress for now)
    if ismember('x',output)
        occ_p=(energyall_p(:,1)<=mu_v);
        bcmap_p=bcmap_p(:);
        if isnan(bcmap_m)
            re.chern_cond=sum(bcmap_p(occ_p))/(2*pi);
        else
            occ_m=(energyall_m(:,1)<=mu_v);
            bcmap_m=bcmap_m(:);
            re.chern_cond=(sum(bcmap_p(occ_p))+sum(bcmap_m(occ_m)))/(2*pi);
        end
        chern_cond_str=sprintf('  \\sigma_{xy}=%.4f',re.chern_cond);
    else
        re.chern_cond=nan;
        chern_cond_str='';
    end
    %% Band structure
    if ismember('f',output)
        re.fig_band=figure;
        
        hold on;
        [energyall_p,energyall_m,~,~]=energyMF_bc(ave1,ave2,'line',epoch,params);
        segment=sqrt(diff(params.k_line(:,1)).^2+diff(params.k_line(:,2)).^2);
        % segment=sqrt(diff(params.k_line_C3(:,1)).^2+diff(params.k_line_C3(:,2)).^2);
        klist=[0;cumsum(segment)];
        for i=1:min(20,size(energyall_p,2))
            scatter(klist,1e3*energyall_p(:,i),5,'r','filled');
            if energyall_m~=0
                scatter(klist,1e3*energyall_m(:,i),5,'b','filled');
            end
        end
        yline(1000*mu_v,'color','k','LineStyle','--');
        yline(1000*mu_c,'color','k','LineStyle','--');
        xline(klist(40),'color','k','LineStyle','--','HandleVisibility','off');
        xline(klist(40+40-1),'color','k','LineStyle','--','HandleVisibility','off');
        xline(klist(40+40-1+80-1),'color','k','LineStyle','--','HandleVisibility','off');
        xticks(klist([1,40,40+40-1,40+40-1+80-1]))
        xticklabels({'\kappa_t','m','\kappa_b','\gamma'});
        xlim([klist(1),klist(40+40-1+80-1)])
        if prod(size(energyall_m))>1
            ylim([1000*min(min(energyall_p(:,1)),min(energyall_m(:,1)))-1,1000*max(energyall_p(:,4*Nai))])
        else
            ylim([1000*min(energyall_p(:,1))-1,1000*max(energyall_p(:,4*Nai))])
        end
        xlabel('|k_m|')
        ylabel('E (meV)')
        drawnow;
        energy_str=sprintf('Gap: %e (meV) E: %e (meV)',re.gap*1000,re.tot*1000);
        epoch_str=sprintf('\nepoch=%d',epoch);
        title(strcat(energy_str,chern_str,chern_cond_str,epoch_str));
        chern_p=re.chern_p;
        chern_m=re.chern_m;
        save(sprintf('energy_w%d_Vb%d_Vz%.1f_epoch%d.mat',params.w*1e3,params.V_b*1e3,params.Vz_t*1e3,epoch),'energyall_p','energyall_m','mu_v','mu_c','klist','chern_p','chern_m')
    else
        re.fig_band=0;
    end
    %% band structure in the BZ
    if ismember('z', output)
        [energyall_p,energyall_m,~,~]=energyMF_bc(ave1,ave2,'bz',epoch,params);
%         figure;scatter(params.k_bz(:,1),params.k_bz(:,2),energyall_m(:,1)*0+5,energyall_m(:,1));
        kx_list=linspace(min(params.k_bz(:,1)),max(params.k_bz(:,1)),1000);
        ky_list=linspace(min(params.k_bz(:,2)),max(params.k_bz(:,2)),1000);
        [kxq,kyq]=meshgrid(kx_list,ky_list);
        energymesh_p=griddata(params.k_bz(:,1),params.k_bz(:,2),energyall_p(:,1),kxq,kyq);
        energymesh_m=griddata(params.k_bz(:,1),params.k_bz(:,2),energyall_m(:,1),kxq,kyq);
        [dos_list,en_list,nu_list,E_vanHove]=DOS(energymesh_p(:),params);

        save(sprintf('DOS_w%d_Vb%d_Vz%.1f_epoch%d.mat',params.w*1e3,params.V_b*1e3,params.Vz_t*1e3,epoch),'energymesh_p','energymesh_m','dos_list','en_list','nu_list','E_vanHove')

    end
    %% density of states
    if ismember('d',output)
        [energyall_p,energyall_m,~,~]=energyMF_bc(ave1,ave2,'bz',epoch,params);
        [dos_list,en_list,nu_list,E_vanHove]=DOS(energyall_p(:,1),params);
    end
    %% Spin
    if ismember('s',output)
        [X_index,Y_index]=meshgrid(linspace(0,1,22));
        [sitesX_index,sitesY_index]=meshgrid(-1:3);
        rmap=[X_index(:),Y_index(:)]*[params.am1;params.am2];
        rsite=[sitesX_index(:),sitesY_index(:)]*[params.aM1;params.aM2];
        re.rmap_x=rmap(:,1);
        re.rmap_y=rmap(:,2);
        re.fig_spin=figure('Position', [10 10 900 600]);
        re.s0_list=zeros(size(rmap));
        re.sx_list=zeros(size(rmap));
        re.sy_list=zeros(size(rmap));
        re.sz_list=zeros(size(rmap));
        for l=1:2
            [s0,sx,sy,sz]=S_r(ave2_n,re.rmap_x,re.rmap_y,l,params);
            re.s0_list(:,l)=s0;
            re.sx_list(:,l)=sx;
            re.sy_list(:,l)=sy;
            re.sz_list(:,l)=sz;
            subplot(2, 2, 2*l-1);
            hold on;
            scatter(re.rmap_x(:)/params.aM,re.rmap_y(:)/params.aM,10,sz(:),'filled');
            quiver(re.rmap_x(:)/params.aM,re.rmap_y(:)/params.aM,(sx(:)/20),(sy(:)/20),'AutoScale','off');
            titlespin=sprintf('%s: epoch=%d\n',(l==1)*'b'+(l==2)*'t',epoch);
            if Nai==3
                ABC_list=[1,(size(rmap,1)-1)/3+1,(size(rmap,1)-1)/3*2+1];
                titlespint='';
                for i=ABC_list
                    % site_label='A'-1+i;
                    titlespint=strcat(titlespint,sprintf('\n(n,\\phi,\\theta)=(%.4f,%.4f,%.4f)',s0(i),atan2(sx(i),sy(i))*180/pi,atan2(sqrt(sx(i)^2+sy(i)^2),sz(i))*180/pi));
                end
            elseif Nai==1
                titlespint=sprintf('A: (n,\\phi,\\theta)=(%.4f,%.4f,%.4f)',s0(1),atan2(sx(1),sy(1))*180/pi,atan2(sqrt(sx(1)^2+sy(1)^2),sz(1))*180/pi);
            end
            title(strcat(titlespin,titlespint));
            daspect([1,1,1]);
            scatter(rsite(:,1)/params.aM,rsite(:,2)/params.aM);
            xlim([min(re.rmap_x)*1.1,max(re.rmap_x)*1.1]/params.aM);
            ylim([min(re.rmap_y)*1.1,max(re.rmap_y)*1.1]/params.aM);
            xlabel('x/|a_M|')
            ylabel('y/|a_M|')
            cb=colorbar;
            cb.Title.String='S_z';
            subplot(2, 2, 2*l);
            hold on;
            scatter(re.rmap_x(:)/params.aM,re.rmap_y(:)/params.aM,10,s0(:),'filled');
            title(sprintf('%s: epoch=%d',(l==1)*'b'+(l==2)*'t',epoch));
            daspect([1,1,1]);
            xlim([min(re.rmap_x)*1.1,max(re.rmap_x)*1.1]/params.aM);
            ylim([min(re.rmap_y)*1.1,max(re.rmap_y)*1.1]/params.aM);
            xlabel('x/|a_M|')
            ylabel('y/|a_M|')
            cb=colorbar;
            cb.Title.String='n';
        end
        
        drawnow
    else
        re.fig_spin=0;
        re.s0_list=nan;
        re.sx_list=nan;
        re.sy_list=nan;
        re.sz_list=nan;
        re.rmap_x=nan;
        re.rmap_y=nan;
    end
    
end