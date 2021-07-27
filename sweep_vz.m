function sweep_vz(nu,epsilon)
    vz_t_list=-38:14;
    chern_p_list=vz_t_list;
    chern_m_list=vz_t_list;
    gap_list=vz_t_list;
    final_list=vz_t_list;
    innergap_list=vz_t_list;
    if nu==[2,1]
        tsymm=1;
    else
        tsymm=0;
    end
    for vz_t_index=1:length(vz_t_list)
        vz_t=vz_t_list(vz_t_index);
        disp(vz_t);
        params=mainTMD('Nmax',3,'V_t',0,'psi_t',240,'V_b',15,'psi_b',-14,'vz_t',vz_t,'vz_b',0,'w',20,'nu',nu,'n',15,'epsilon',epsilon,'tsymm',tsymm);
        [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(0,0,params);
        if nu==[2,1]
            [ave1,ave2,occ]=average(energyall,wfall,0,1,params); %for QAHE
        else
            [ave1,ave2,occ]=average(energyall,wfall,1,0,params); %for QSHE
        end

        Nk=size(params.k,1);
        tot=totalenergy(V1_ave_delta,V2_ave_delta,ave1,ave2,params);
        tot_list=[tot];
        
        for i=1:100
            [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,params);    
            [ave1_n,ave2_n,occ]=average(energyall,wfall,0,1,params);
            tot=totalenergy(V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,params);
            tot_list(end+1)=tot;
            ave1=ave1_n;
            ave2=ave2_n;
            if abs(tot_list(end)-tot_list(end-1))<1e-7
                break;
            end
        end

        energyall_sort=sort(energyall(:));
        gap=energyall_sort(Nk*params.nu(1)/(params.nu(2))+1)-energyall_sort(Nk*params.nu(1)/(params.nu(2)));
        innergap=energyall_sort(Nk*params.nu(1)/(params.nu(2)))-energyall_sort(Nk*params.nu(1)/(params.nu(2))-1);
        finali=i;
        [chern_p,chern_m]=chern_gs(ave1,ave2,2,params);

        chern_p_list(vz_t_index)=chern_p(1);
        chern_m_list(vz_t_index)=chern_m(1);
        gap_list(vz_t_index)=gap;
        final_list(vz_t_index)=tot;
        innergap_list(vz_t_index)=innergap;
        finali_list(vz_t_index)=finali;
    end

    save(sprintf('phase_nu%d,%d_ep%d.mat',params.nu(1),params.nu(2),epsilon),'chern_p_list','chern_m_list','gap_list','final_list','innergap_list','finali_list')
    
    