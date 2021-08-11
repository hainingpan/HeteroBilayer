function tot=totalenergy(V1_ave_delta,V2_ave_delta,ave1_n,ave2_n,epoch,params)
b_set=(params.b);
q_set=(params.q);
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
T_tensor=reshape(T,[Nq,Nb,2,2,Nq,Nb,2,2,Nk]); %q_b,b_b,l_b,t_b,q_a,b_a,l_a,t_a,k_b
%ave2_n k_b,q_b,b_b,l_b,t_b,q_a,b_a,l_a,t_a
% H0=sum(permute(ave2_n,[2,3,4,5,6,7,8,9,1]).*T_tensor,'all');
H0=collapse(permute(ave2_n,[2,3,4,5,6,7,8,9,1]).*T_tensor,1:length(size(T_tensor)));

A=Nk*Nai*params.area;
if isa(V1_ave_delta,'double') && V1_ave_delta==0
    H1=0;
else
    % V1=params.V1;   % q_g,q_d,b_g,b_d
    %ave1_a: q_g,q_d,b_g,b_d
    % V1_ave=V1.*ave1; %q_g,q_d,b_g,b_d
    % delta=params.delta_tensor1; %q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d
    % V1_ave_delta=ttt2(V1_ave,delta,[1,2,3,4],[3,4,7,8],[],[]);  %q_a,q_b,b_a,b_b
    %ave1_b: q_b,q_a,b_b,b_a
    H1=collapse(permute(ave1_n,[2,1,4,3]).*V1_ave_delta,1:length(size(V1_ave_delta)))/(2*A*params.epsilon);
end
if isa(V2_ave_delta,'double') && V2_ave_delta==0
    H2=0;
else
    % V2=params.V2;  % k_a,k_b,q_a,q_d,b_a,b_d 
    %ave2_a: k_a,q_d,b_d,l_a,t_a,q_g,b_g,l_b,t_b
    % V2_ave=ttt2(V2,ave2,[1],[1],[4,6],[2,3]);   %: q_d,b_d,k_b,q_a,b_a,l_a,t_a,q_g,b_g,l_b,t_b
    % delta=params.delta_tensor2; %q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d
    % V2_ave_delta=ttt2(V2_ave,delta,[1,2,8,9],[4,8,3,7],[4,5],[1,5]); % q_a,b_a,k_b,l_a,t_a,l_b,t_b,q_b,b_b
    %ave2_b: k_b,q_b,b_b,l_b,t_b,q_a,b_a,l_a,t_a
    H2=collapse(permute(V2_ave_delta,[3,8,9,6,7,1,2,4,5]).*ave2_n,1:length(size(ave2_n)))/(2*A*params.epsilon);
end

tot=H0+H1-H2;
tot_err=abs(imag(tot));
assert(tot_err<1e-12,sprintf("hermitian error: %e\n",tot_err));
tot=real(tot)/(Nk*Nai); %energy per unit cell

