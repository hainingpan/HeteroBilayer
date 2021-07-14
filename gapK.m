Eglist=linspace(-50,-30,40);
gapKlist=Eglist*0;
parfor index=1:length(Eglist)
    Eg=Eglist(index);
    params_mK=mainTMD('valley',-1,'vz_t',Eg);
    [en,vec]=energyTMD(params_mK.kt(1),params_mK.kt(2),params_mK);
    gapKlist(index)=en(1)-en(2);
end

m=(params_mK.kb+params_mK.kt)/2;
kt_m_x=linspace(params_mK.kt(1),m(1),40);
kt_m_y=linspace(params_mK.kt(2),m(2),40);
m_kb_x=linspace(m(1),params_mK.kb(1),40);
m_kb_y=linspace(m(2),params_mK.kb(2),40);
kb_gamma_x=linspace(params_mK.kb(1),0,40);
kb_gamma_y=linspace(params_mK.kb(2),0,40);


kxlist=[kt_m_x,m_kb_x,kb_gamma_x];
kylist=[kt_m_y,m_kb_y,kb_gamma_y];

gapminlist=Eglist*0;
parfor index=1:length(Eglist)
    Eg=Eglist(index);
    params_mK=mainTMD('valley',-1,'vz_t',Eg);
    enlist=zeros(length(kxlist),1);
    for index2=1:length(kxlist)
        [en,~]=energyTMD(kxlist(index2),kylist(index2),params_mK);
        enlist(index2)=en(1)-en(2);
    end
    gaplist(index)=min(enlist);
end

chern1list=Eglist*0;
chern2list=Eglist*0;
parfor index=1:length(Eglist)
    Eg=Eglist(index);
    params_mK=mainTMD('valley',-1,'vz_t',Eg);
    [~,~,~,~,~,~,chern1]=berrycurvature(1,params_mK);
    [~,~,~,~,~,~,chern2]=berrycurvature(2,params_mK);
    chern1list(index)=chern1;
    chern2list(index)=chern2;
end





