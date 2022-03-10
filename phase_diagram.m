% function phase_diagram()
% 
% vz_list=linspace(-50,50,101);
% w_list=linspace(0,30,61);

vz_list=linspace(-50,0,101);
w_list=linspace(0,20,81);
topo_pK_map=zeros(length(vz_list),length(w_list));
gap_indirect_map=zeros(length(vz_list),length(w_list));
gap_direct_map=zeros(length(vz_list),length(w_list));
Vb=7;
for vz_i=1:length(vz_list)
    vz=vz_list(vz_i);
    for w_i=1:length(w_list)
        w=w_list(w_i);
        params_pK=mainTMD('Nmax',2,'V_t',0,'psi_t',240,'V_b',Vb,'psi_b',-14,'vz_t',vz,'vz_b',0,'w',w);
        [kcxmap,kcymap,kcx2map,kcy2map,bcmap,omega,chern_pK]=berrycurvature(1,params_pK);
        topo_pK_map(vz_i,w_i)=chern_pK;
        [gap_indirect_map(vz_i,w_i),gap_direct_map(vz_i,w_i)]=find_band_gap(params_pK);
    end
end
figure;surf(w_list,vz_list,topo_pK_map,'edgecolor','none');view(2);xlabel('w (meV)');ylabel('V_{zt} (meV)');colorbar;

% save('phase_diagram.mat','gap_indirect_map','gap_direct_map','topo_pK_map','vz_list','w_list')
save(sprintf('pd_single_particle_Vb%d.mat',Vb),'gap_indirect_map','gap_direct_map','topo_pK_map','vz_list','w_list')




        
