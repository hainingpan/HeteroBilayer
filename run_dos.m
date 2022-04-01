vz_list=-(0:5:50);
for vz=vz_list
    disp(vz)
    dos(vz)
end


function dos(vz)
    
Nmax=2;
nu=[1,1];
w=12;
vz=vz;
Vb=7;
ep=12;

params=mainTMD('SDW',10e-3,'d',5,'Nmax',Nmax,'V_t',0,'psi_t',240,'V_b',Vb,'psi_b',-14,'vz_t',vz,'vz_b',0,'w',w,'nu',nu,'n',21,'epsilon',ep,'shift',1);
[energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(0,0,0,params);
[ave1,ave2,occ]=average(energyall,wfall,0,params); 

[gap,tot]=plotline_2(energyall,0,0,V1_ave_delta,V2_ave_delta,ave1,ave2,0,strcat('fsd',params.chern),params);
end