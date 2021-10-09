vz_list=linspace(-100,0,51);
topo_pK_list=vz_list;
gap_list=vz_list;
% topo_mK_list=vz_list;
for vz_i=1:length(vz_list)
    vz=vz_list(vz_i);
    params_pK=mainTMD('Nmax',2,'V_t',0,'psi_t',240,'V_b',7,'psi_b',-14,'vz_t',vz,'vz_b',0,'w',12);
    % params_pK=mainTMD('Nmax',3,'valley',-1,'V_t',0,'psi_t',240,'V_b',4.1,'psi_b',14,'vz_t',vz,'vz_b',0,'w',1.3);
    [kcxmap,kcymap,kcx2map,kcy2map,bcmap,omega,chern_pK]=berrycurvature(1,params_pK);
%     [kcxmap,kcymap,kcx2map,kcy2map,bcmap,omega,chern_mK]=berrycurvature(1,params_mK);
    topo_pK_list(vz_i)=chern_pK;
%     topo_mK_list(vz_i)=chern_mK;
    gap_list(vz_i)=find_band_gap(params_pK);
end

% figure;plot(vz_list,1000*max(0,gap_list));
figure;plot(vz_list,1000*gap_list);
hold on;
plot(vz_list,topo_pK_list);