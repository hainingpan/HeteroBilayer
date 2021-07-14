psi_list=linspace(180,300,30);
vz_t_list=linspace(-30,0,10);
chern1_list=zeros(length(psi_list),length(vz_t_list));
chern2_list=zeros(length(psi_list),length(vz_t_list));
chern3_list=zeros(length(psi_list),length(vz_t_list));

for index=1:length(psi_list)
    for i2=1:length(vz_t_list)
        params_mK=mainTMD('valley',-1,'V_t',4.1,'V_b',4.1,'psi_t',psi_list(index),'psi_b',14,'vz_t',vz_t_list(i2),'vz_b',0,'w',6);
        [kcxmap,kcymap,kcx2map,kcy2map,bcmap,omega,chern1]=berrycurvature(1,params_mK);
        [kcxmap,kcymap,kcx2map,kcy2map,bcmap,omega,chern2]=berrycurvature(2,params_mK);
        [kcxmap,kcymap,kcx2map,kcy2map,bcmap,omega,chern3]=berrycurvature(3,params_mK);
        chern1_list(index,i2)=chern1;
        chern2_list(index,i2)=chern2;
        chern3_list(index,i2)=chern3;
    end
end

% figure;plot(psi_list,chern1_list);
% hold on;
% plot(psi_list,chern2_list);
% plot(psi_list,chern3_list);