function [energyall_p,energyall_m,wfall_p,wfall_m]=energyMF_bc(ave1,ave2,tag,params)
b_set=(params.b);
q_set=(params.q);
Nb=size(b_set,1);
Nq=size(q_set,1);
if strcmp(tag,'bc')
    k_beta_set=params.k_bc;
end
if strcmp(tag,'line')
    k_beta_set=params.k_line;
    % k_beta_set=[0,0];
end
Nk0=size(params.k,1);
Nk=size(k_beta_set,1);
kb=params.kb;
kt=params.kt;
m_b=params.m_b;
m_t=params.m_t;
Vz_b=params.Vz_b;
Vz_t=params.Vz_t;
%H0
q_mat=eye(Nq);
b_mat=eye(Nb);
energyall_p=zeros(Nk,Nb*Nq*2);
energyall_m=zeros(Nk,Nb*Nq*2);
% valley_index=zeros(Nk,Nb*Nq*2*2);
wfall_p=zeros(Nk,Nb*Nq*2,Nb*Nq*2);
wfall_m=zeros(Nk,Nb*Nq*2,Nb*Nq*2);
T=zeros(Nb*Nq*2*2,Nb*Nq*2*2,Nk);

alpha=0.00729735; % eV*nm
[k_a_x,k_b_x,q_a_x,q_d_x,b_a_x,b_d_x]=ndgrid(params.k(:,1),k_beta_set(:,1),params.q(:,1),params.q(:,1),params.b(:,1),params.b(:,1));
[k_a_y,k_b_y,q_a_y,q_d_y,b_a_y,b_d_y]=ndgrid(params.k(:,2),k_beta_set(:,2),params.q(:,2),params.q(:,2),params.b(:,2),params.b(:,2));
q_abs=sqrt((k_a_x+b_d_x+q_d_x-k_b_x-b_a_x-q_a_x).^2+(k_a_y+b_d_y+q_d_y-k_b_y-b_a_y-q_a_y).^2);
qd=q_abs*params.d;
V2=alpha*2*pi*tanh(qd+1e-18)./(qd+1e-18)*params.d;


for k_beta_index=1:Nk
    kx=k_beta_set(k_beta_index,1);
    ky=k_beta_set(k_beta_index,2);
    [qx_set_grid,bx_set_grid]=ndgrid(q_set(:,1),b_set(:,1));
    [qy_set_grid,by_set_grid]=ndgrid(q_set(:,2),b_set(:,2));

    k_b_p_x=kx+qx_set_grid(:)+bx_set_grid(:)-kb(1);
    k_b_p_y=ky+qy_set_grid(:)+by_set_grid(:)-kb(2);

    k_t_p_x=kx+qx_set_grid(:)+bx_set_grid(:)-kt(1);
    k_t_p_y=ky+qy_set_grid(:)+by_set_grid(:)-kt(2);

    k_b_m_x=-kx-qx_set_grid(:)-bx_set_grid(:)-kb(1);
    k_b_m_y=-ky-qy_set_grid(:)-by_set_grid(:)-kb(2);

    k_t_m_x=-kx-qx_set_grid(:)-bx_set_grid(:)-kt(1);
    k_t_m_y=-ky-qy_set_grid(:)-by_set_grid(:)-kt(2);

    T_b_m=-1/(2*m_b)*diag(k_b_m_x.^2+k_b_m_y.^2);
    T_t_m=-1/(2*m_t)*diag(k_t_m_x.^2+k_t_m_y.^2);
    T_m=[T_b_m,0*T_b_m;0*T_b_m,T_t_m];

    T_b_p=-1/(2*m_b)*diag(k_b_p_x.^2+k_b_p_y.^2);
    T_t_p=-1/(2*m_t)*diag(k_t_p_x.^2+k_t_p_y.^2);
    T_p=[T_b_p,0*T_b_p;0*T_b_p,T_t_p];
    
    T_tau=[T_p,0*T_m;0*T_m,T_m];

    Delta_b_p=kron(params.Delta_b_p,q_mat);
    Delta_t_p=kron(params.Delta_t_p,q_mat);
    Delta_T_p=kron(params.Delta_T_p,q_mat);
    Delta_TT_p=kron(params.Delta_TT_p,q_mat);
    Delta_b_m=kron(params.Delta_b_m,q_mat);
    Delta_t_m=kron(params.Delta_t_m,q_mat);
    Delta_T_m=kron(params.Delta_T_m,q_mat);
    Delta_TT_m=kron(params.Delta_TT_m,q_mat);

    Delta_p=[Delta_b_p,Delta_T_p;Delta_TT_p,Delta_t_p];
    Delta_m=[Delta_b_m,Delta_T_m;Delta_TT_m,Delta_t_m];

    Delta_tau=[Delta_p,0*Delta_p;0*Delta_p,Delta_m];

    Vz_b_mat=Vz_b*kron(b_mat,q_mat);
    Vz_t_mat=Vz_t*kron(b_mat,q_mat);
    Vz_mat=[Vz_b_mat,0*Vz_b_mat;0*Vz_b_mat,Vz_t_mat];
    Vz_tau=[Vz_mat,0*Vz_mat;0*Vz_mat,Vz_mat];
    H0=T_tau+Delta_tau+Vz_tau;
    T(:,:,k_beta_index)=-H0.';
end

A=Nk0*Nq*params.area;
if ave1==0
    H1=0;
else
    V1=params.V1;   % q_g,q_d,b_g,b_d
    %ave1: q_g,q_d,b_g,b_d
    V1_ave=V1.*ave1; %q_g,q_d,b_g,b_d
    delta=params.delta_tensor1; %q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d
    V1_ave_delta=ttt2(V1_ave,delta,[1,2,3,4],[3,4,7,8],[],[]);  %q_a,q_b,b_a,b_b
    hartree=reshape(permute(V1_ave_delta,[1,3,2,4]),Nq*Nb,Nq*Nb); %q_a,b_a,q_b,b_b
    H1=(kron(eye(4),hartree));
    H1=repmat(H1,[1,1,Nk])/(A*params.epsilon);
end

if ave2==0
    H2=0;
else
    % V2=params.V2;  % k_a,k_b,q_a,q_d,b_a,b_d 
    %ave2: k_a,q_d,b_d,l_a,t_a,q_g,b_g,l_b,t_b
    V2_ave=ttt2(V2,ave2,[1],[1],[4,6],[2,3]);   %: q_d,b_d,k_b,q_a,b_a,l_a,t_a,q_g,b_g,l_b,t_b
    delta=params.delta_tensor2; %q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d
    V2_ave_delta=ttt2(V2_ave,delta,[1,2,8,9],[4,8,3,7],[4,5],[1,5]); % q_a,b_a,k_b,l_a,t_a,l_b,t_b,q_b,b_b
    H2=permute(V2_ave_delta,[1,2,4,5,8,9,6,7,3]); %q_a,b_a,l_a,t_a,q_b,b_b,l_b,t_b,k_b
    H2=reshape(H2,[Nq*Nb*2*2,Nq*Nb*2*2,Nk])/(A*params.epsilon);
end
H=T+H1-H2;

herr=max(sum(abs(H-conj(permute(H,[2,1,3]))),[1,2]));
assert(herr<1e-12,sprintf("hermitian error: %e\n",herr));
H=1/2*(H+conj(permute(H,[2,1,3])));
for k_beta_index=1:Nk
    H_tau=H(:,:,k_beta_index);
    if isequal(H_tau(1:end/2,end/2+1:end),0*H_tau(1:end/2,end/2+1:end))
        H_p=H_tau(1:end/2,1:end/2);
        H_m=H_tau(end/2+1:end,end/2+1:end);
        [vec_p,val_p]=eig(H_p);
        [vec_m,val_m]=eig(H_m);
        val_p=real(diag(val_p));
        val_m=real(diag(val_m));
        [val_p,I_p]=sort(val_p);
        [val_m,I_m]=sort(val_m);
        vec_p=vec_p(:,I_p);
        vec_m=vec_m(:,I_m);
    else
        [vec,val]=eig(H_tau);
        val=real(diag(val));
        [val,I]=sort(val);
        vec=vec(:,I);
        valley_index=0;
    end
    energyall_p(k_beta_index,:)=val_p;
    energyall_m(k_beta_index,:)=val_m;
    % valley_index(k_beta_index,:)=valley;

    for ii=1:Nb*Nq*2
        wfall_p(k_beta_index,ii,:)=vec_p(:,ii); %k,n,component for +K
        wfall_m(k_beta_index,ii,:)=vec_m(:,ii); %k,n,component for -K
    end
end
end
