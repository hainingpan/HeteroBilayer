function [gap,tot]=plotline_2(energyall,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,params)

    energyall_sort=sort(energyall(:));
    Nk=size(params.k,1);
    Nq=size(params.q,1);
    mu=energyall_sort(Nk*Nq*params.nu(1)/(params.nu(2)));
    mu1=energyall_sort(1+Nk*Nq*params.nu(1)/(params.nu(2)));
    
    gap=energyall_sort(Nk*Nq*params.nu(1)/(params.nu(2))+1)-energyall_sort(Nk*Nq*params.nu(1)/(params.nu(2)));
    tot=totalenergy(V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,params);

    kxlist=params.k_line(:,1);
    kylist=params.k_line(:,2);
    segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
    klist=[0;cumsum(segment)];
    
    figure;
    chern_p=[0];
    chern_m=[0];
    % [chern_p,chern_m]=chern_gs(ave1,ave2,1,epoch,params);

    title(sprintf('Gap: %e (meV) E: %e (meV)\n +K:{%.4f}\n-K:{%.4f}',gap*1000,tot*1000,chern_p,chern_m));
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
    ylim([1000*min(energyall_p(:,1)),1000*max(energyall_p(:,4*Nq))])
    drawnow;
    