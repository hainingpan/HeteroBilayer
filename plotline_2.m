function [gap,tot,fig_band,fig_spin,chern_p,chern_m]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,output,params)

    Nai=size(params.ailist,1);  % The expansion of super cell
    if ismember('g',output)
        [energyall_p,energyall_m,~,~]=energyMF_bc(ave1,ave2,'dense',epoch,params);
        if prod(size(energyall_m))>1
            energyall_sort=sort([energyall_p(:);energyall_m(:)]);
        else
            energyall_sort=sort(energyall_p(:));
        end
        Nk=size(params.k_dense,1);
    else
        Nk=size(params.k,1);
        energyall_sort=sort(energyall(:));
    end
    mu_v=energyall_sort(Nk*Nai*params.nu(1)/(params.nu(2)));
    mu_c=energyall_sort(1+Nk*Nai*params.nu(1)/(params.nu(2)));
    % mu_v_index=find(energyall==mu_v,1);
    % mu_c_index=find(energyall==mu_c,1);
    % [mu_v_k_index,mu_v_l_index]=ind2sub(size(energyall),mu_v_index);
    % [mu_c_k_index,mu_c_l_index]=ind2sub(size(energyall),mu_c_index);
    % gap_str=sprintf('\nv:%d,%d c:%d,%d',mu_v_k_index,mu_v_l_index,mu_c_k_index,mu_c_l_index);
    gap=mu_c-mu_v;
    % tot=totalenergy(V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,params);
    if prod(size(V1_ave_delta))>1 && prod(size(V2_ave_delta))>1
        tot=totalenergy(V1_ave_delta,V2_ave_delta,ave1,ave2,epoch,params);
    else
        tot=totalenergy(V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,params);
    end

    %% Chern number
    if ismember('c',output)
        [chern_p,chern_m]=chern_gs(ave1,ave2,epoch,params);
        chern_str=sprintf('\n+K:{%.4f} -K:{%.4f}',chern_p,chern_m);
    else
        chern_p=nan;
        chern_m=nan;
        chern_str='';
    end
    
    
    %% Band structure
    if ismember('f',output)
        fig_band=figure;
        energy_str=sprintf('Gap: %e (meV) E: %e (meV)',gap*1000,tot*1000);
        epoch_str=sprintf('\nepoch=%d',epoch);
        title(strcat(energy_str,chern_str,epoch_str));
        hold on;
        [energyall_p,energyall_m,~,~]=energyMF_bc(ave1,ave2,'line',epoch,params);
        segment=sqrt(diff(params.k_line(:,1)).^2+diff(params.k_line(:,2)).^2);
        klist=[0;cumsum(segment)];
        for i=1:20
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
        ylim([1000*min(energyall_p(:,1))-1,1000*max(energyall_p(:,4*Nai))])
        xlabel('|k_m|')
        ylabel('E (meV)')
        drawnow;
    else
        fig_band=0;
    end
    %% Spin
    if ismember('s',output)
        [X_index,Y_index]=meshgrid(linspace(0,1,22));
        [sitesX_index,sitesY_index]=meshgrid(-1:3);
        rmap=[X_index(:),Y_index(:)]*[params.am1;params.am2];
        rsite=[sitesX_index(:),sitesY_index(:)]*[params.aM1;params.aM2];
        rmap_x=rmap(:,1);
        rmap_y=rmap(:,2);
        [s0,sx,sy,sz]=S_r(ave2_n,rmap_x,rmap_y,1,params);
        fig_spin=figure;
        hold on;
        scatter(rmap_x(:)/params.aM,rmap_y(:)/params.aM,10,sz(:),'filled')
        quiver(rmap_x(:)/params.aM,rmap_y(:)/params.aM,(sx(:)),(sy(:)));
        title(sprintf('epoch=%d',epoch));
        daspect([1,1,1]);
        scatter(rsite(:,1)/params.aM,rsite(:,2)/params.aM);
        xlim([min(rmap_x)*1.1,max(rmap_x)*1.1]/params.aM);
        ylim([min(rmap_y)*1.1,max(rmap_y)*1.1]/params.aM);
        xlabel('x/|a_M|')
        ylabel('y/|a_M|')
        cb=colorbar;
        cb.Title.String='S_z';
        drawnow
    else
        fig_spin=0;
    end
    
end