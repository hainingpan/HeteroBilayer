figure;
% +K
% w=12;
% psi_b=-14;
% V_b=7;
% vz_t=-20;

w=0;
vz_t=0;
V_b=6;
psi_b=0;

%  params_pK=mainTMD('valley',1,'vz_t',0);
% params_pK=mainTMD('Nmax',2,'valley',1,'V_t',4.1,'psi_t',240,'V_b',4.1,'psi_b',-14,'vz_t',-30,'vz_b',0,'w',1.3);
params_pK=mainTMD('Nmax',3,'valley',1,'V_b',V_b,'psi_b',psi_b,'vz_t',vz_t,'vz_b',0,'w',w);

% params_pK=mainTMD('Nmax',3,'valley',1,'d',100e-9*5.076e6,'a_t',3.42e-10*5.076e6,'a_b',3.575e-10*5.076e6,'V_t',20.1,'psi_t',240,'V_b',20.1,'psi_b',-14,'vz_t',-100,'vz_b',0,'w',0,'nu',[1,1],'n',15,'epsilon',20);

kxlist=params_pK.k_line(:,1);
kylist=params_pK.k_line(:,2);
segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
klist=[0;cumsum(segment)];

energy_p_list=zeros(2*length(params_pK.h1index),length(kxlist)); %initialize
for i=1:length(kxlist)
    energy_p_list(:,i)=energyTMD(kxlist(i),kylist(i),params_pK);
end
plot(klist,1000*energy_p_list,'r')
hold on

%% -K
%  params_mK=mainTMD('valley',-1,'vz_t',-40);
params_mK=mainTMD('Nmax',3,'valley',-1,'V_b',V_b,'psi_b',psi_b,'vz_t',vz_t,'vz_b',0,'w',w);
% params_mK=mainTMD('Nmax',3,'m_b',1.2,'valley',-1,'V_t',0,'psi_t',240,'V_b',15,'psi_b',-14,'vz_t',0,'vz_b',0,'w',20);
% params_mK=mainTMD('Nmax',3,'valley',-1,'d',100e-9*5.076e6,'a_t',3.42e-10*5.076e6,'a_b',3.575e-10*5.076e6,'V_t',20.1,'psi_t',240,'V_b',20.1,'psi_b',-14,'vz_t',-100,'vz_b',0,'w',0,'nu',[1,1],'n',15,'epsilon',20);

 %   params_mK=mainTMD('valley',-1,'a_t',3.28e-10*5.076e6,'a_b',3.28e-10*5.076e6,'m_b',0.45,'m_t',0.45,'V_b',20,'V_t',20,'psi_b',-108,'psi_t',108,'w',20,'Vz_b',0,'Vz_t',0,'omega',1,'theta',3);

energy_m_list=zeros(2*length(params_mK.h1index),length(kxlist)); %initialize
for i=1:length(kxlist)
    energy_m_list(:,i)=energyTMD(kxlist(i),kylist(i),params_mK);
end
plot(klist,1000*energy_m_list,'b')

xline(klist(40),'color','k','LineStyle','--','HandleVisibility','off');
xline(klist(40+40-1),'color','k','LineStyle','--','HandleVisibility','off');
xline(klist(40+40-1+80-1),'color','k','LineStyle','--','HandleVisibility','off');


xticks(klist([1,40,40+40-1,40+40-1+80-1]))
xticklabels({'\kappa_t','m','\kappa_b','\gamma'});
xlim([klist(1),klist(end)])
ylim([min(1000*energy_m_list(10,:)),1.2*max(1000*energy_m_list(1,:))])
% ylim([0,1.2*max(1000*energylist(1,:))])

save('bandstructure.mat','energy_p_list','energy_m_list','w','V_b','psi_b','vz_t','klist','kxlist','kylist','segment');