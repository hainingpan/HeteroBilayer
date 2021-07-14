function [val,vec,H]=energyTMD(kx,ky,params)
%energy at (kx,ky)
if params.valley==1
    DeltaTmat=params.DeltaTmat;
    DeltaTTmat=params.DeltaTTmat;
    Deltabmat=params.Deltabmat;
    Deltatmat=params.Deltatmat;
    Vz_b=params.Vz_b;
    Vz_t=params.Vz_t;
else
    DeltaTmat=conj(params.DeltaTmat);
    DeltaTTmat=conj(params.DeltaTTmat);
    Deltabmat=conj(params.Deltabmat);
    Deltatmat=conj(params.Deltatmat);
    Vz_b=conj(params.Vz_b);
    Vz_t=conj(params.Vz_t);
end
m_b=params.m_b;
m_t=params.m_t;
kb=params.kb;
kt=params.kt;
bM1=params.bM1;
bM2=params.bM2;
h1index=params.h1index;
h2index=params.h2index;
kblist=(params.valley)*[kx,ky]-kb+h1index(:)*bM1+h2index(:)*bM2;
ktlist=(params.valley)*[kx,ky]-kt+h1index(:)*bM1+h2index(:)*bM2;
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