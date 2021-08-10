function tot=totalenergy_2(V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,params)
% Gauge transformation to shift the momentum to bottom -> k^2
b_set=(params.b);
q_set=(params.q);
Nb=size(b_set,1);
Nq=size(q_set,1);
k_beta_set=params.k;
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
energyall=zeros(Nk,Nb*Nq*2*2);
valley_index=zeros(Nk,Nb*Nq*2*2);
wfall=zeros(Nk,Nb*Nq*2*2,Nb*Nq*2*2);
T=zeros(Nb*Nq*2*2,Nb*Nq*2*2,Nk);
for k_beta_index=1:Nk
    kx=k_beta_set(k_beta_index,1);
    ky=k_beta_set(k_beta_index,2);
    [qx_set_grid,bx_set_grid]=ndgrid(q_set(:,1),b_set(:,1));
    [qy_set_grid,by_set_grid]=ndgrid(q_set(:,2),b_set(:,2));

    k_b_p_x=kx+qx_set_grid(:)+bx_set_grid(:);
    k_b_p_y=ky+qy_set_grid(:)+by_set_grid(:);

    k_b_m_x=kx+qx_set_grid(:)+bx_set_grid(:);
    k_b_m_y=ky+qy_set_grid(:)+by_set_grid(:);

    k_t_p_x=kx+qx_set_grid(:)+bx_set_grid(:)-(kt(1)-kb(1));
    k_t_p_y=ky+qy_set_grid(:)+by_set_grid(:)-(kt(2)-kb(2));

    k_t_m_x=kx+qx_set_grid(:)+bx_set_grid(:)+(kt(1)-kb(1));
    k_t_m_y=ky+qy_set_grid(:)+by_set_grid(:)+(kt(2)-kb(2));

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
T_tensor=reshape(T,[Nq,Nb,2,2,Nq,Nb,2,2,Nk]); %q_b,b_b,l_b,t_b,q_a,b_a,l_a,t_a,k_b
%ave2_n k_b,q_b,b_b,l_b,t_b,q_a,b_a,l_a,t_a
H0=sum(permute(ave2_n,[2,3,4,5,6,7,8,9,1]).*T_tensor,'all');

A=Nk*Nq*params.area;
if V1_ave_delta==0
    H1=0;
else
    % V1=params.V1;   % q_g,q_d,b_g,b_d
    %ave1_a: q_g,q_d,b_g,b_d
    % V1_ave=V1.*ave1; %q_g,q_d,b_g,b_d
    % delta=params.delta_tensor1; %q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d
    % V1_ave_delta=ttt2(V1_ave,delta,[1,2,3,4],[3,4,7,8],[],[]);  %q_a,q_b,b_a,b_b
    %ave1_b: q_b,q_a,b_b,b_a
    H1=sum(permute(ave1_n,[2,1,4,3]).*V1_ave_delta,'all')/(2*A*params.epsilon);
end
if V2_ave_delta==0
    H2=0;
else
    % V2=params.V2;  % k_a,k_b,q_a,q_d,b_a,b_d 
    %ave2_a: k_a,q_d,b_d,l_a,t_a,q_g,b_g,l_b,t_b
    % V2_ave=ttt2(V2,ave2,[1],[1],[4,6],[2,3]);   %: q_d,b_d,k_b,q_a,b_a,l_a,t_a,q_g,b_g,l_b,t_b
    % delta=params.delta_tensor2; %q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d
    % V2_ave_delta=ttt2(V2_ave,delta,[1,2,8,9],[4,8,3,7],[4,5],[1,5]); % q_a,b_a,k_b,l_a,t_a,l_b,t_b,q_b,b_b
    %ave2_b: k_b,q_b,b_b,l_b,t_b,q_a,b_a,l_a,t_a
    H2=sum(permute(V2_ave_delta,[3,8,9,6,7,1,2,4,5]).*ave2_n,'all')/(2*A*params.epsilon);
end

tot=H0+H1-H2;
tot_err=abs(imag(tot));
assert(tot_err<1e-12,sprintf("hermitian error: %e\n",tot_err));
tot=real(tot)/(Nk*Nb); %energy per unit cell

