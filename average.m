function [ave1,ave2,occ]=average(energyall,wfall,epoch,params)
    k_beta_set=(params.k);
    Nk=size(k_beta_set,1);
    q_set=(params.q);
    Nq=size(q_set,1);
    Nai=size(params.ailist,1);
    b_set=params.b;
    Nb=size(b_set,1);

    energyall_sort=sort(energyall(:));
    if epoch==0
        if params.fermisurface==1
            mu=energyall_sort(Nk*Nai*params.nu(1)/(params.nu(2)));
            occ=(energyall<=mu); %k,n
        else
            occ=energyall*0;
            occ(:,1:2)=1;
            occ(:,3:end)=0;
        end
        v_p=sum(abs(wfall(:,:,1:end/2)).^2,3);
        v_m=sum(abs(wfall(:,:,end/2+1:end)).^2,3);
        if params.valley_polarized==1
            occ=(v_p>v_m).*occ;
        elseif params.valley_polarized==-1        
            occ=~((v_p>v_m)).*occ;
        end
    else
        mu=energyall_sort(Nk*Nai*params.nu(1)/(params.nu(2)));
        occ=(energyall<=mu); %k,n
    end

    d=reshape(wfall,Nk,Nq*Nb*4,Nq,Nb,2,2); % d_{k,n,q,b,l,tau}
    d_conj=conj(d);
    occ_d=repmat(occ,[1,1,Nq,Nb,2,2]).*d;   %k_a,n,q_d,b_d,l_a,t_a
    ave1=ttt2(d_conj,occ_d,[1,2,5,6],[1,2,5,6],[],[]); %q_g,b_g,q_d,b_d

    ave1=permute(ave1,[1,3,2,4]); %q_g,q_d,b_g,b_d

    %check hermittian?

    ave2=ttt2(d_conj,occ_d,[2],[2],[1],[1]); %k_a,q_d,b_d,l_a,t_a,q_g,b_g,l_b,t_b
    if params.tsymm==1
        ave2_2=conj(ave2(end:-1:1,end:-1:1,end:-1:1,:,end:-1:1,end:-1:1,end:-1:1,:,end:-1:1));
        ave2=(ave2+ave2_2)/2;
        % ave2_2=conj(ave2(end:-1:1,end:-1:1,end:-1:1,:,end:-1:1,end:-1:1,end:-1:1,:,end:-1:1));
        % ave2_err=max(ave2-ave2_2,[],'all');
    end
    occ=double(occ);