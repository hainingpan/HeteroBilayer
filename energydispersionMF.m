figure;
 params_mK=mainTMD('Nmax',3,'valley',1,'V_t',4.1,'psi_t',240,'V_b',4.1,'psi_b',-14,'vz_t',-30,'vz_b',0,'w',1.3,'nu',[2,1]);

[energyall,wfall]=energyMF(0,0,params_mK);
klist=params_mK.klist;
energylist=energyall';
plot(klist,1000*energylist,'k')

line(klist([40,40]),[-10,200],'color','k','LineStyle','--','HandleVisibility','off');
line(klist([80,80]),[-10,200],'color','k','LineStyle','--','HandleVisibility','off');
line(klist([120,120]),[-10,200],'color','k','LineStyle','--','HandleVisibility','off');


xticks(klist([1,40,80,120]))
xticklabels({'\kappa_t','m','\kappa_b','\gamma'});
xlim([klist(1),klist(end)])
% ylim([min(1000*energylist(end-4,:)),1.2*max(1000*energylist(end,:))])
% ylim([1.2*max(1000*energylist(1,:)),])
ylim([-10,130])