% params=mainTMD('Nmax',3,'d',100e-9*5.076e6,'a_t',3.42e-10*5.076e6,'a_b',3.575e-10*5.076e6,'V_t',20.1,'psi_t',240,'V_b',20.1,'psi_b',-14,'vz_t',-30,'vz_b',0,'w',1.3,'nu',[1,1],'n',15,'epsilon',5);
params=mainTMD('Nmax',2,'d',5e-9*5.076e6,'V_t',4.1,'psi_t',240,'V_b',4.1,'psi_b',-14,'vz_t',-30,'vz_b',0,'w',1.3,'nu',[2,1],'n',21,'epsilon',2);
[energyall,wfall,valley_index]=energyMF(0,0,params);
[ave1,ave2,occ,~]=average(energyall,wfall,0,1,params);

gap=plotline(energyall,occ,valley_index,params);
gap_list=[gap];

fig1=figure;
for i=1:100
    [energyall,wfall,valley_index]=energyMF(ave1,ave2,params);
    [ave1,ave2,occ,~]=average(energyall,wfall,0,1,params);
    gap=plotline(energyall,occ,valley_index,params);
    gap_list(end+1)=gap;
    figure(fig1);
    plot(gap_list*1e3);
end

