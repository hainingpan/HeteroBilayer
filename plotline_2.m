function [gap,tot]=plotline_2(energyall,wfall,occ,valley_index,ave1,ave2,V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,params)

    energyall_sort=sort(energyall(:));
    Nk=size(params.k,1);
    mu=energyall_sort(Nk*params.nu(1)/(params.nu(2)));
    mu1=energyall_sort(1+Nk*params.nu(1)/(params.nu(2)));
    
    gap=energyall_sort(Nk*params.nu(1)/(params.nu(2))+1)-energyall_sort(Nk*params.nu(1)/(params.nu(2)));
    tot=totalenergy(V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,params);
    % kt=params.bM1*1/3+params.bM2*2/3;
    % kb=params.bM1*2/3+params.bM2*1/3;
    % m=(kb+kt)/2;
    % gamma=[0,0];
    % kt_m_x=linspace(kt(1),m(1),40);
    % kt_m_y=linspace(kt(2),m(2),40);
    % m_kb_x=linspace(m(1),kb(1),40);
    % m_kb_y=linspace(m(2),kb(2),40);
    % kb_gamma_x=linspace(kb(1),gamma(1),40);
    % kb_gamma_y=linspace(kb(2),gamma(2),40);
    % kxlist=[kt_m_x(1:end-1),m_kb_x(1:end-1),kb_gamma_x];
    % kylist=[kt_m_y(1:end-1),m_kb_y(1:end-1),kb_gamma_y];

    kxlist=params.k_line(:,1);
    kylist=params.k_line(:,2);
    segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
    klist=[0;cumsum(segment)];
    
    figure;
    [chern_p,chern_m]=chern_gs(ave1,ave2,2,params);
    title(sprintf('Gap: %e (meV) E: %e (meV)\n +K:{%.4f,%.4f}\n-K:{%.4f,%.4f}',gap*1000,tot*1000,chern_p,chern_m));
    hold on;
    [energyall_p,energyall_m,~,~]=energyMF_bc(ave1,ave2,'line',params);
    for i=1:10
        scatter(klist,1e3*energyall_p(:,i),5,'r','filled');
        scatter(klist,1e3*energyall_m(:,i),5,'b','filled');
    end

    % for i=1:10
    %     en_q=griddata(params.k(:,1),params.k(:,2),1000*energyall(:,i),kxlist,kylist);
    %     valley_q=griddata(params.k(:,1),params.k(:,2),valley_index(:,i),kxlist,kylist);
    %     occ_q=5*griddata(params.k(:,1),params.k(:,2),occ(:,i),kxlist,kylist)+1;
    %     valley_q_rgb=zeros(length(kxlist),3);
    %     for j=1:length(kxlist)
    %         valley_q_rgb(j,:)=[(valley_q(j)+1)/2,0,(-valley_q(j)+1)/2];
    %     end
    %     scatter(klist,en_q,occ_q*10,valley_q_rgb,'filled');
    %     % plot(klist,en_q,'b')
    % end
    
    yline(1000*mu,'color','k','LineStyle','--');
    yline(1000*mu1,'color','k','LineStyle','--');
    xline(klist(40),'color','k','LineStyle','--','HandleVisibility','off');
    xline(klist(79),'color','k','LineStyle','--','HandleVisibility','off');
    xline(klist(118),'color','k','LineStyle','--','HandleVisibility','off');
    % line(klist([1,118]),[1000*mu,1000*mu],'color','k','LineStyle','--');
    % line(klist([40,40]),[-400,200],'color','k','LineStyle','--','HandleVisibility','off');
    % line(klist([79,79]),[-400,200],'color','k','LineStyle','--','HandleVisibility','off');
    % line(klist([118,118]),[-200,200],'color','k','LineStyle','--','HandleVisibility','off');
    
    
    xticks(klist([1,40,79,118]))
    xticklabels({'\kappa_t','m','\kappa_b','\gamma'});
    xlim([klist(1),klist(118)])
    ylim([1000*min(energyall_p(:,1)),1000*max(energyall_p(:,4))])
    
    drawnow;
    