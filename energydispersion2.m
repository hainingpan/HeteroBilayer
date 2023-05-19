parameters=mainTMD('valley',-1,'a_t',3.28e-10*5.076e6,'a_b',3.28e-10*5.076e6,'m_b',0.45,'m_t',0.45,'V_b',20,'V_t',20,'psi_b',108,'psi_t',-108,'w',20,'Vz_b',0,'Vz_t',0,'omega',1,'theta',3);

parameters.kpp=[-parameters.kb(1),parameters.kb(2)];

kpp_gamma_x=linspace(parameters.kpp(1),0,40);
kpp_gamma_y=linspace(parameters.kpp(2),0,40);
gamma_kn_x=linspace(0,parameters.kt(1),40);
gamma_kn_y=linspace(0,parameters.kt(2),40);
kn_kp_x=linspace(parameters.kt(1),parameters.kb(1),40);
kn_kp_y=linspace(parameters.kt(2),parameters.kb(2),40);
kp_kpp_x=linspace(parameters.kb(1),parameters.kpp(1),80);
kp_kpp_y=linspace(parameters.kb(2),parameters.kpp(2),80);



kxlist=[kpp_gamma_x,gamma_kn_x,kn_kp_x,kp_kpp_x];
kylist=[kpp_gamma_y,gamma_kn_y,kn_kp_y,kp_kpp_y];

segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
klist=[0,cumsum(segment)];

energylist=zeros(2*(2*parameters.Nmax+1)^2,length(kxlist)); %initialize
for i=1:length(kxlist)
    energylist(:,i)=energyTMD(kxlist(i),kylist(i),parameters);
end
figure;
plot(klist,1000*energylist)
hold on
line(klist([40,40]),[-400,30],'color','k','LineStyle','--','HandleVisibility','off');
line(klist([80,80]),[-400,30],'color','k','LineStyle','--','HandleVisibility','off');
line(klist([120,120]),[-200,30],'color','k','LineStyle','--','HandleVisibility','off');


xticks(klist([1,40,80,120,200]))
xticklabels({'\kappa_+^\prime','\gamma','\kappa_-','\kappa_+','\kappa_+^\prime'});
xlim([klist(1),klist(end)])
ylim([min(1000*energylist(5,:)),1.2*max(1000*energylist(1,:))])