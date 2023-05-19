function [energyall,wfall,valley_index,V1_ave_delta,V2_ave_delta]=energyMF(ave1,ave2,epoch,params)
b_set=params.b;
q_set=params.q;
Nb=size(b_set,1);
Nq=size(q_set,1);
Nai=size(params.ailist,1);  % The expansion of super cell
k_beta_set=params.k;
Nk=size(k_beta_set,1);

Vz_b=params.Vz_b;
Vz_t=params.Vz_t;
%H0
q_mat=eye(Nq);
b_mat=eye(Nb);
energyall=zeros(Nk,Nb*Nq*2*2);
valley_index=zeros(Nk,Nb*Nq*2*2);
wfall=zeros(Nk,Nb*Nq*2*2,Nb*Nq*2*2);

if epoch>0
    Delta_b_p=(params.Delta_b_p);
    Delta_t_p=(params.Delta_t_p);
    Delta_T_p=(params.Delta_T_p);
    Delta_TT_p=(params.Delta_TT_p);
    Delta_b_m=(params.Delta_b_m);
    Delta_t_m=(params.Delta_t_m);
    Delta_T_m=(params.Delta_T_m);
    Delta_TT_m=(params.Delta_TT_m);
else
    Delta_b_p=(params.Delta_b_0_p);
    Delta_t_p=(params.Delta_t_0_p);
    Delta_T_p=(params.Delta_T_0_p);
    Delta_TT_p=(params.Delta_TT_0_p);
    Delta_b_m=(params.Delta_b_0_m);
    Delta_t_m=(params.Delta_t_0_m);
    Delta_T_m=(params.Delta_T_0_m);
    Delta_TT_m=(params.Delta_TT_0_m);
end

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

A=Nk*Nai*params.area;
if isa(ave1,'double') && ave1==0
    V1_ave_delta=0;
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
    V2_ave_delta=0;
    H2=0;
else
    V2=params.V2;  % k_a,k_b,q_a,q_d,b_a,b_d 
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
clear V1_ave V2_ave  delta hartree

herr=max(sum(abs(H-conj(permute(H,[2,1,3]))),[1,2]));
assert(herr<1e-12,sprintf("hermitian error: %e\n",herr));
H=1/2*(H+conj(permute(H,[2,1,3])));
% check_TR(T,params);
% if numel(ave1)>1 && numel(ave2)>1
%     check_TR(H1,params);
%     check_TR(H2,params);
% end
for k_beta_index=1:Nk
    H_tau=H(:,:,k_beta_index);
    if isequal(H_tau(1:end/2,end/2+1:end),0*H_tau(1:end/2,end/2+1:end))
        H_p=H_tau(1:end/2,1:end/2);
        H_m=H_tau(end/2+1:end,end/2+1:end);
        [vec_p,val_p]=eig(H_p);
        [vec_m,val_m]=eig(H_m);
        val_p=real(diag(val_p));
        val_m=real(diag(val_m));
        vec_p_1=[vec_p;0*vec_p];
        vec_m_1=[0*vec_m;vec_m];
        val=[val_p;val_m];
        valley=[ones(size(val_p));-ones(size(val_m))];
        vec=[vec_p_1,vec_m_1];
        [val,I]=sort(val);
        vec=vec(:,I);
        valley=valley(I);
        valley_index(k_beta_index,:)=valley;
    else
        [vec,val]=eig(H_tau);
        val=real(diag(val));
        [val,I]=sort(val);
        vec=vec(:,I);
        valley_index=0;
    end
    energyall(k_beta_index,:)=val;
    
    for ii=1:Nb*Nq*2*2
        wfall(k_beta_index,ii,:)=vec(:,ii); %k,n,component
    end
end

end
