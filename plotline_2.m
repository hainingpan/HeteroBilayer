function [gap,tot,fig_band,fig_spin]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,output,params)

    energyall_sort=sort(energyall(:));
    Nk=size(params.k,1);
    % Nq=size(params.q,1);
    Nai=size(params.ailist,1);  % The expansion of super cell
    mu=energyall_sort(Nk*Nai*params.nu(1)/(params.nu(2)));
    mu1=energyall_sort(1+Nk*Nai*params.nu(1)/(params.nu(2)));
    
    gap=energyall_sort(Nk*Nai*params.nu(1)/(params.nu(2))+1)-energyall_sort(Nk*Nai*params.nu(1)/(params.nu(2)));
    tot=totalenergy(V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,params);

    segment=sqrt(diff(params.k_line(:,1)).^2+diff(params.k_line(:,2)).^2);
    klist=[0;cumsum(segment)];
    %% Chern number

    if ismember('c',output)
        [chern_p,chern_m]=chern_gs(ave1,ave2,1,epoch,params);
        chern_str=sprintf('\n+K:{%.4f}\n-K:{%.4f}',chern_p,chern_m);
    else
        chern_str='';
    end

    %% Band structure
    if ismember('f',output)
        fig_band=figure;
        energy_str=sprintf('Gap: %e (meV) E: %e (meV)',gap*1000,tot*1000);
        title(strcat(energy_str,chern_str));
        hold on;
        [energyall_p,energyall_m,~,~]=energyMF_bc(ave1,ave2,'line',epoch,params);
        for i=1:10
            scatter(klist,1e3*energyall_p(:,i),5,'r','filled');
            if energyall_m~=0
                scatter(klist,1e3*energyall_m(:,i),5,'b','filled');
            end
        end
        yline(1000*mu,'color','k','LineStyle','--');
        yline(1000*mu1,'color','k','LineStyle','--');
        xline(klist(40),'color','k','LineStyle','--','HandleVisibility','off');
        xline(klist(40+40-1),'color','k','LineStyle','--','HandleVisibility','off');
        xline(klist(40+40-1+80-1),'color','k','LineStyle','--','HandleVisibility','off');
        xticks(klist([1,40,40+40-1,40+40-1+80-1]))
        xticklabels({'\kappa_t','m','\kappa_b','\gamma'});
        xlim([klist(1),klist(40+40-1+80-1)])
        ylim([1000*min(energyall_p(:,1)),1000*max(energyall_p(:,4*Nai))])
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