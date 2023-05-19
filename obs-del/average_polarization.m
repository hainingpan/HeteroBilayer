function [ave1,ave2,occ,valley_index]=average_polarization(energyall,wfall,valley,params)
    % obsolete
    k_beta_set=(params.k);
    Nk=size(k_beta_set,1);
    energyall_sort=sort(energyall(:));
    mu=energyall_sort(Nk*2);
    occ=(energyall<mu); %k,n

    v_p=sum(abs(wfall(:,:,1:end/2)).^2,3);
    v_m=sum(abs(wfall(:,:,end/2+1:end)).^2,3);
    valley_index=(v_p-v_m); % 1: +K valley; -1: -K valley
    if valley==1
        occ=(v_p>v_m).*occ;
    else
        occ=~((v_p>v_m)).*occ;
    end

    q_set=(params.q);
    Nq=size(q_set,1);
    b_set=params.b;
    Nb=size(b_set,1);
    % N=size(energyall,1);

    d=reshape(wfall,Nk,Nq*Nb*4,Nq,Nb,2,2); % d_{k,n,q,b,l,tau}
    d_conj=conj(d);
    occ_d=repmat(occ,[1,1,Nq,Nb,2,2]).*d;
    ave1=ttt2(d_conj,occ_d,[1,2],[1,2],[5,6],[5,6]); %l_a,t_a,q_g,b_g,q_d,b_d

    ave1=permute(ave1,[3,4,5,6,1,2]); %q_g,b_g,q_d,b_d,l_a,t_a

    %check hermittian?

    ave2=ttt2(d_conj,occ_d,[2],[2],[1],[1]); %k_a,q_d,b_d,l_a,t_a,q_g,b_g,l_b,t_b

    occ=double(occ);









