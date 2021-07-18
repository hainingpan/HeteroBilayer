figure;
% +K
w=1.3;
%  params_pK=mainTMD('valley',1,'vz_t',0);
% params_pK=mainTMD('Nmax',2,'valley',1,'V_t',4.1,'psi_t',240,'V_b',4.1,'psi_b',-14,'vz_t',-30,'vz_b',0,'w',1.3);
params_pK=mainTMD('Nmax',2,'d',5e-9*5.076e6,'valley',1,'V_t',4.1,'psi_t',240,'V_b',4.1,'psi_b',-14,'vz_t',-30,'vz_b',0,'w',1.3,'nu',[2,1],'n',21);

% params_pK=mainTMD('valley',1,'w',0,'V_t',4.1,'psi_t',240);

m=(params_pK.kb+params_pK.kt)/2;
% gamma=-params_pK.shift;
gamma=[0,0];
kt_m_x=linspace(params_pK.kt(1),m(1),40);
kt_m_y=linspace(params_pK.kt(2),m(2),40);
m_kb_x=linspace(m(1),params_pK.kb(1),40);
m_kb_y=linspace(m(2),params_pK.kb(2),40);
kb_gamma_x=linspace(params_pK.kb(1),gamma(1),40);
kb_gamma_y=linspace(params_pK.kb(2),gamma(2),40);


kxlist=[kt_m_x,m_kb_x,kb_gamma_x];
kylist=[kt_m_y,m_kb_y,kb_gamma_y];

segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
klist=[0,cumsum(segment)];

energylist=zeros(2*length(params_pK.h1index),length(kxlist)); %initialize
for i=1:length(kxlist)
    energylist(:,i)=energyTMD(kxlist(i),kylist(i),params_pK);
end
plot(klist,1000*energylist,'b')
hold on

%% -K
%  params_mK=mainTMD('valley',-1,'vz_t',-40);
 params_mK=mainTMD('Nmax',2,'valley',-1,'V_t',4.1,'psi_t',240,'V_b',4.1,'psi_b',-14,'vz_t',-30,'vz_b',0,'w',1.3);
%   params_mK=mainTMD('valley',-1,'a_t',3.28e-10*5.076e6,'a_b',3.28e-10*5.076e6,'m_b',0.45,'m_t',0.45,'V_b',20,'V_t',20,'psi_b',-108,'psi_t',108,'w',20,'Vz_b',0,'Vz_t',0,'omega',1,'theta',3);

m=(params_mK.kb+params_mK.kt)/2;
% gamma=-params_mK.kb0;
gamma=[0,0];
kt_m_x=linspace(params_mK.kt(1),m(1),40);
kt_m_y=linspace(params_mK.kt(2),m(2),40);
m_kb_x=linspace(m(1),params_mK.kb(1),40);
m_kb_y=linspace(m(2),params_mK.kb(2),40);
kb_gamma_x=linspace(params_mK.kb(1),gamma(1),40);
kb_gamma_y=linspace(params_mK.kb(2),gamma(2),40);


kxlist=[kt_m_x,m_kb_x,kb_gamma_x];
kylist=[kt_m_y,m_kb_y,kb_gamma_y];

segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
klist=[0,cumsum(segment)];

energylist=zeros(2*length(params_mK.h1index),length(kxlist)); %initialize
for i=1:length(kxlist)
    energylist(:,i)=energyTMD(kxlist(i),kylist(i),params_mK);
end
plot(klist,1000*energylist,'r')


line(klist([40,40]),[-400,30],'color','k','LineStyle','--','HandleVisibility','off');
line(klist([80,80]),[-400,30],'color','k','LineStyle','--','HandleVisibility','off');
line(klist([120,120]),[-200,30],'color','k','LineStyle','--','HandleVisibility','off');


xticks(klist([1,40,80,120]))
xticklabels({'\kappa_t','m','\kappa_b','\gamma'});
xlim([klist(1),klist(end)])
% ylim([min(1000*energylist(3,:)),1.2*max(1000*energylist(1,:))])
ylim([-130,1.2*max(1000*energylist(1,:))])