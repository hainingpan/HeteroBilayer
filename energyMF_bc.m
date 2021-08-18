function [energyall_p,energyall_m,wfall_p,wfall_m]=energyMF_bc(ave1,ave2,tag,epoch,params)
b_set=(params.b);
q_set=(params.q);
Nb=size(b_set,1);
Nq=size(q_set,1);
Nai=size(params.ailist,1);  % The expansion of super cell
if strcmp(tag,'bc')
    k_beta_set=params.k_bc;
end
if strcmp(tag,'line')
    k_beta_set=params.k_line;
    % k_beta_set=params.k_line_C3;
    % k_beta_set=[0,0];
end
if strcmp(tag,'dense')
    k_beta_set=params.k_dense;
end
Nk0=size(params.k,1);
Nk=size(k_beta_set,1);
Vz_b=params.Vz_b;
Vz_t=params.Vz_t;
%H0
q_mat=eye(Nq);
b_mat=eye(Nb);


alpha=0.00729735; % eV*nm
[k_a_x,k_b_x,q_a_x,q_d_x,b_a_x,b_d_x]=ndgrid(params.k(:,1),k_beta_set(:,1),params.q(:,1),params.q(:,1),params.b(:,1),params.b(:,1));
k_x=k_a_x+b_d_x+q_d_x-k_b_x-b_a_x-q_a_x;
clear k_a_x b_d_x q_d_x k_b_x b_a_x q_a_x 
[k_a_y,k_b_y,q_a_y,q_d_y,b_a_y,b_d_y]=ndgrid(params.k(:,2),k_beta_set(:,2),params.q(:,2),params.q(:,2),params.b(:,2),params.b(:,2));
k_y=k_a_y+b_d_y+q_d_y-k_b_y-b_a_y-q_a_y;
clear k_a_y b_d_y q_d_y k_b_y b_a_y q_a_y
qd=params.d*sqrt(k_x.^2+k_y.^2);
clear k_x k_y
V2=alpha*2*pi*tanh(qd+1e-18)./(qd+1e-18)*params.d;
clear qd;
V2=tensor(V2,[Nk0,Nk,Nq,Nq,Nb,Nb]);
Delta_b_p=(params.Delta_b_p);
Delta_t_p=(params.Delta_t_p);
Delta_T_p=(params.Delta_T_p);
Delta_TT_p=(params.Delta_TT_p);
Delta_b_m=(params.Delta_b_m);
Delta_t_m=(params.Delta_t_m);
Delta_T_m=(params.Delta_T_m);
Delta_TT_m=(params.Delta_TT_m);

Delta_p=[Delta_b_p,Delta_T_p;Delta_TT_p,Delta_t_p];
Delta_m=[Delta_b_m,Delta_T_m;Delta_TT_m,Delta_t_m];

Delta_tau=[Delta_p,0*Delta_p;0*Delta_p,Delta_m];

Vz_b_mat=Vz_b*kron(b_mat,q_mat);
Vz_t_mat=Vz_t*kron(b_mat,q_mat);
Vz_mat=[Vz_b_mat,0*Vz_b_mat;0*Vz_b_mat,Vz_t_mat];
Vz_tau=[Vz_mat,0*Vz_mat;0*Vz_mat,Vz_mat];
if epoch==0
    S_tau=params.S_tau;
else
    S_tau=0;
end
shift=params.shift;
T0=T_gen(k_beta_set,shift,params);
T=T0+repmat(Delta_tau,[1,1,Nk])+repmat(Vz_tau,[1,1,Nk])+repmat(S_tau,[1,1,Nk]);

T=-permute(T,[2,1,3]);

A=Nk0*Nai*params.area;
if isa(ave1,'double') && ave1==0
    H1=0;
else
    V1=params.V1;   % q_g,q_d,b_g,b_d
    %ave1: q_g,q_d,b_g,b_d
    V1_ave=V1.*ave1; %q_g,q_d,b_g,b_d
    delta=params.delta_tensor1; %q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d
    V1_ave_delta=ttt2(V1_ave,delta,[1,2,3,4],[3,4,7,8],[],[]);  %q_a,q_b,b_a,b_b
    hartree=reshape(permute(V1_ave_delta,[2,4,1,3]),[Nq*Nb,Nq*Nb]); %{q_b,b_b},{q_a,b_a}
    % hartree=reshape(permute(V1_ave_delta,[1,3,2,4]),[Nq*Nb,Nq*Nb]); %{q_b,b_b},{q_a,b_a}
    H1=(kron(eye(4),hartree.data));
    H1=repmat(H1,[1,1,Nk])/(A*params.epsilon);
    % H1=permute(H1,[2,1,3]);
end

if isa(ave2,'double') && ave2==0
    H2=0;
else
    % V2=params.V2;  % k_a,k_b,q_a,q_d,b_a,b_d 
    %ave2: k_a,q_d,b_d,l_a,t_a,q_g,b_g,l_b,t_b
    V2_ave=ttt2(V2,ave2,[1],[1],[4,6],[2,3]);   %: q_d,b_d,k_b,q_a,b_a,l_a,t_a,q_g,b_g,l_b,t_b
    delta=params.delta_tensor2; %q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d
    V2_ave_delta=ttt2(V2_ave,delta,[1,2,8,9],[4,8,3,7],[4,5],[1,5]); % q_a,b_a,k_b,l_a,t_a,l_b,t_b,q_b,b_b
    H2=permute(V2_ave_delta,[8,9,6,7,1,2,4,5,3]); %q_b,b_b,l_b,t_b,q_a,b_a,l_a,t_a,k_b
    % H2=permute(V2_ave_delta,[1,2,4,5,8,9,6,7,3]); %q_a,b_a,l_a,t_a,q_b,b_b,l_b,t_b,k_b
    H2=reshape(H2,[Nq*Nb*2*2,Nq*Nb*2*2,Nk])/(A*params.epsilon);
    H2=H2.data;
    % H2=permute(H2,[2,1,3]);
end
H=T+H1-H2;
% clear V1_ave V2_ave  delta hartree

herr=max(sum(abs(H-conj(permute(H,[2,1,3]))),[1,2]));
assert(herr<1e-12,sprintf("hermitian error: %e\n",herr));
H=1/2*(H+conj(permute(H,[2,1,3])));

if isequal(H(1:end/2,end/2+1:end,1),0*H(1:end/2,end/2+1:end,1))
    energyall_p=zeros(Nk,Nb*Nq*2);
    energyall_m=zeros(Nk,Nb*Nq*2);
    wfall_p=zeros(Nk,Nb*Nq*2,Nb*Nq*2);
    wfall_m=zeros(Nk,Nb*Nq*2,Nb*Nq*2);
else
    energyall_p=zeros(Nk,Nb*Nq*4);
    energyall_m=0;
    wfall_p=zeros(Nk,Nb*Nq*4,Nb*Nq*4);
    wfall_m=0;
end
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
        energyall_p(k_beta_index,:)=val_p;
        energyall_m(k_beta_index,:)=val_m;
        for ii=1:Nb*Nq*2
            wfall_p(k_beta_index,ii,:)=vec_p(:,ii); %k,n,component for +K
            wfall_m(k_beta_index,ii,:)=vec_m(:,ii); %k,n,component for -K
        end
    else
        [vec,val]=eig(H_tau);
        val=real(diag(val));
        [val,I]=sort(val);
        vec=vec(:,I);
        energyall_p(k_beta_index,:)=val;
        for ii=1:Nb*Nq*2
            wfall_p(k_beta_index,ii,:)=vec(:,ii); %k,n,component for both valleys
        end
    end

end
end
