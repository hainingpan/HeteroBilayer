function [energyall,wfall]=energyMF(ave,Vave,params)
Nb=size(params.b,1);
Nq=size(params.q,1);
k_beta_set=(params.k);
N=size(k_beta_set,1);
b_set=(params.b);
q_set=(params.q);
kb=params.kb;
kt=params.kt;
m_b=params.m_b;
m_t=params.m_t;
Vz_b=params.Vz_b;
Vz_t=params.Vz_t;
%H0
q_mat=eye(Nq);
b_mat=eye(Nb);
energyall=zeros(N,Nb*Nq*2*2);
wfall=zeros(N,Nb*Nq*2*2,Nb*Nq*2*2);
T=zeros(Nb*Nq*2*2,Nb*Nq*2*2,N);
for k_beta_index=1:N
    kx=k_beta_set(k_beta_index,1);
    ky=k_beta_set(k_beta_index,2);
    [qx_set_grid,bx_set_grid]=ndgrid(q_set(:,1),b_set(:,1));
    [qy_set_grid,by_set_grid]=ndgrid(q_set(:,2),b_set(:,2));

    k_b_m_x=kx+qx_set_grid(:)+bx_set_grid(:)-kb(1);
    k_b_m_y=ky+qy_set_grid(:)+by_set_grid(:)-kb(2);

    k_t_m_x=kx+qx_set_grid(:)+bx_set_grid(:)-kt(1);
    k_t_m_y=ky+qy_set_grid(:)+by_set_grid(:)-kt(2);

    k_b_p_x=-kx+qx_set_grid(:)+bx_set_grid(:)-kb(1);
    k_b_p_y=-ky+qy_set_grid(:)+by_set_grid(:)-kb(2);

    k_t_p_x=-kx+qx_set_grid(:)+bx_set_grid(:)-kt(1);
    k_t_p_y=-ky+qy_set_grid(:)+by_set_grid(:)-kt(2);

    T_b_m=-1/(2*m_b)*diag(k_b_m_x.^2+k_b_m_y.^2);
    T_t_m=-1/(2*m_t)*diag(k_t_m_x.^2+k_t_m_y.^2);
    T_m=[T_b_m,0*T_b_m;0*T_b_m,T_t_m];

    T_b_p=-1/(2*m_b)*diag(k_b_p_x.^2+k_b_p_y.^2);
    T_t_p=-1/(2*m_t)*diag(k_t_p_x.^2+k_t_p_y.^2);
    T_p=[T_b_p,0*T_b_p;0*T_b_p,T_t_p];
    
    T_tau=[T_m,0*T_m;0*T_m,T_p];

    Deltabmat=kron(params.Deltabmat,q_mat);
    Deltatmat=kron(params.Deltatmat,q_mat);
    DeltaTmat=kron(params.DeltaTmat,q_mat);
    DeltaTTmat=kron(params.DeltaTTmat,q_mat);
    Delta=[Deltabmat,DeltaTmat;DeltaTTmat,Deltatmat];
    Delta_tau=[Delta,0*Delta;0*Delta,conj(Delta)];

    Vz_b_mat=Vz_b*kron(b_mat,q_mat);
    Vz_t_mat=Vz_t*kron(b_mat,q_mat);
    Vz_mat=[Vz_b_mat,0*Vz_b_mat;0*Vz_b_mat,Vz_t_mat];
    Vz_tau=[Vz_mat,0*Vz_mat;0*Vz_mat,conj(Vz_mat)];
    H0=T_tau+Delta_tau+Vz_tau;
    T(:,:,k_beta_index)=-H0.';
end

if ave==0
    H1=0;
else

end

if Vave==0
    H2=0;
else
    
end
H=T-H1+H2;



herr=max(sum(abs(H-conj(permute(H,[2,1,3]))),[1,2]));
assert(herr<1e-12,sprintf("hermitian error: %e\n",herr));
H=1/2*(H+conj(permute(H,[2,1,3])));
for k_beta_index=1:N
   [vec,val]=eig(H(:,:,k_beta_index));
    val=real(diag(val));
    [val,I]=sort(val);
    vec=vec(:,I);
    energyall(k_beta_index,:)=val;
    for ii=1:Nb*Nq*2*2
        wfall(k_beta_index,ii,:)=vec(:,ii);
    end
end
