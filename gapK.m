Eglist=linspace(-50,-30,40);
gapKlist=Eglist*0;
parfor index=1:length(Eglist)
    Eg=Eglist(index);
    params_mK=mainTMD('valley',-1,'vz_t',Eg);
    [en,vec]=energyTMD(params_mK.kt(1),params_mK.kt(2),params_mK);
    gapKlist(index)=en(1)-en(2);
end

