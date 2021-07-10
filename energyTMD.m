function [val,vec]=energyTMD(kx,ky,parameters)
%energy at (kx,ky)
if parameters.valley==-1
    DeltaTmat=parameters.DeltaTmat;
    DeltaTTmat=parameters.DeltaTTmat;
    Deltabmat=parameters.Deltabmat;
    Deltatmat=parameters.Deltatmat;
    Vz_b=parameters.Vz_b;
    Vz_t=parameters.Vz_t;
else
    DeltaTmat=conj(parameters.DeltaTmat);
    DeltaTTmat=conj(parameters.DeltaTTmat);
    Deltabmat=conj(parameters.Deltabmat);
    Deltatmat=conj(parameters.Deltatmat);
    Vz_b=conj(parameters.Vz_b);
    Vz_t=conj(parameters.Vz_t);
end
m_b=parameters.m_b;
m_t=parameters.m_t;
kb=parameters.kb;
kt=parameters.kt;
bM1=parameters.bM1;
bM2=parameters.bM2;
h1index=parameters.h1index;
h2index=parameters.h2index;
kblist=(-parameters.valley)*[kx,ky]-kb+h1index(:)*bM1+h2index(:)*bM2;
ktlist=(-parameters.valley)*[kx,ky]-kt+h1index(:)*bM1+h2index(:)*bM2;
H11=-1/(2*m_b)*diag(dot(kblist,kblist,2))+Deltabmat+Vz_b*eye(length(Deltabmat));
H22=-1/(2*m_t)*diag(dot(ktlist,ktlist,2))+Deltatmat+Vz_t*eye(length(Deltabmat));
H12=DeltaTmat;
H21=DeltaTTmat;
H=[H11,H12;H21,H22];
[vec,val]=eig(H);
val=diag(val);
val=val(end:-1:1);
vec=vec(:,end:-1:1);
end