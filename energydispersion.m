%% +K
%  params_pK=mainTMD('valley',1,'vz_t',0);
% params_pK=mainTMD('valley',1);

% m=(params_pK.kb+params_pK.kt)/2;
% kt_m_x=linspace(params_pK.kt(1),m(1),40);
% kt_m_y=linspace(params_pK.kt(2),m(2),40);
% m_kb_x=linspace(m(1),params_pK.kb(1),40);
% m_kb_y=linspace(m(2),params_pK.kb(2),40);
% kb_gamma_x=linspace(params_pK.kb(1),0,40);
% kb_gamma_y=linspace(params_pK.kb(2),0,40);
% 
% 
% kxlist=[kt_m_x,m_kb_x,kb_gamma_x];
% kylist=[kt_m_y,m_kb_y,kb_gamma_y];
% 
% segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
% klist=[0,cumsum(segment)];

% energylist=zeros(2*(2*params_pK.Nmax+1)^2,length(kxlist)); %initialize
% for i=1:length(kxlist)
%     energylist(:,i)=energyTMD(kxlist(i),kylist(i),params_pK);
% end
figure;
% plot(klist,1000*energylist,'b')
% hold on

%% -K
params_mK=mainTMD('valley',-1,'vz_t',00);
% params_mK=mainTMD('valley',-1);

m=(params_mK.kb+params_mK.kt)/2;
kt_m_x=linspace(params_mK.kt(1),m(1),40);
kt_m_y=linspace(params_mK.kt(2),m(2),40);
m_kb_x=linspace(m(1),params_mK.kb(1),40);
m_kb_y=linspace(m(2),params_mK.kb(2),40);
kb_gamma_x=linspace(params_mK.kb(1),0,40);
kb_gamma_y=linspace(params_mK.kb(2),0,40);


kxlist=[kt_m_x,m_kb_x,kb_gamma_x];
kylist=[kt_m_y,m_kb_y,kb_gamma_y];

segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
klist=[0,cumsum(segment)];

energylist=zeros(2*(2*params_mK.Nmax+1)^2,length(kxlist)); %initialize
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
ylim([min(1000*energylist(3,:)),1.2*max(1000*energylist(1,:))])