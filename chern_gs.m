function [chern_p,chern_m]=chern_gs(ave1,ave2,epoch,params)
% ignore q for now. ignore nonabelian wilson loop
    Nai=size(params.ailist,1);  % The expansion of super cell

    [~,~,wfall_p_0,wfall_m_0]=energyMF_bc(ave1,ave2,'bc',epoch,params);
    if wfall_m_0==0
        n=sqrt(size(wfall_p_0,1));
        l_index=1:Nai;
        wfall_p=wfall_p_0(:,l_index,:);
        umap_p=permute(reshape(wfall_p,[n,n,Nai,size(wfall_p,3)]),[3,4,1,2]);
        nindex=1:n;
        nindex_inc=mod(nindex,n)+1;
        u1_p=umap_p(:,:,nindex,nindex); %Nai,NL,n,n,
        u2_p=umap_p(:,:,nindex_inc,nindex);
        u3_p=umap_p(:,:,nindex_inc,nindex_inc);
        u4_p=umap_p(:,:,nindex,nindex_inc);
        bcmap_p=zeros(n,n);
        for i=nindex
            for j=nindex
                U12=conj(u1_p(:,:,i,j))*u2_p(:,:,i,j).';
                U23=conj(u2_p(:,:,i,j))*u3_p(:,:,i,j).';
                U34=conj(u3_p(:,:,i,j))*u4_p(:,:,i,j).';
                U41=conj(u4_p(:,:,i,j))*u1_p(:,:,i,j).';
                bcmap_p(i,j)=angle(det(U12*U23*U34*U41));
            end
        end
        chern_p=sum(bcmap_p(:))/(2*pi);
        chern_m=nan;
    else
        n=sqrt(size(wfall_p_0,1));
        l_index=1:Nai;
        wfall_p=wfall_p_0(:,l_index,:);
        wfall_m=wfall_m_0(:,l_index,:);
        umap_p=permute(reshape(wfall_p,[n,n,Nai,size(wfall_p,3)]),[3,4,1,2]);
        umap_m=permute(reshape(wfall_m,[n,n,Nai,size(wfall_m,3)]),[3,4,1,2]);

        nindex=1:n;
        nindex_inc=mod(nindex,n)+1;
        u1_p=umap_p(:,:,nindex,nindex); %Nai,NL,n,n,
        u2_p=umap_p(:,:,nindex_inc,nindex);
        u3_p=umap_p(:,:,nindex_inc,nindex_inc);
        u4_p=umap_p(:,:,nindex,nindex_inc);
        bcmap_p=zeros(n,n);
        for i=nindex
            for j=nindex
                U12=conj(u1_p(:,:,i,j))*u2_p(:,:,i,j).';
                U23=conj(u2_p(:,:,i,j))*u3_p(:,:,i,j).';
                U34=conj(u3_p(:,:,i,j))*u4_p(:,:,i,j).';
                U41=conj(u4_p(:,:,i,j))*u1_p(:,:,i,j).';
                bcmap_p(i,j)=angle(det(U12*U23*U34*U41));
            end
        end

        u1_m=umap_m(:,:,nindex,nindex);
        u2_m=umap_m(:,:,nindex_inc,nindex);
        u3_m=umap_m(:,:,nindex_inc,nindex_inc);
        u4_m=umap_m(:,:,nindex,nindex_inc);


        bcmap_m=zeros(n,n);
        for i=nindex
            for j=nindex
                U12=conj(u1_m(:,:,i,j))*u2_m(:,:,i,j).';
                U23=conj(u2_m(:,:,i,j))*u3_m(:,:,i,j).';
                U34=conj(u3_m(:,:,i,j))*u4_m(:,:,i,j).';
                U41=conj(u4_m(:,:,i,j))*u1_m(:,:,i,j).';
                bcmap_m(i,j)=angle(det(U12*U23*U34*U41));
            end
        end

        chern_p=sum(bcmap_p(:))/(2*pi);
        chern_m=sum(bcmap_m(:))/(2*pi);
    end
end


