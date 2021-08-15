function [chern_p,chern_m]=chern_gs(ave1,ave2,epoch,parameters)
        % ignore q for now. ignore nonabelian wilson loop
        Nai=size(params.ailist,1);  % The expansion of super cell

        [~,~,wfall_p_0,wfall_m_0]=energyMF_bc(ave1,ave2,'bc',epoch,parameters);

        if wfall_m_0==0
                plot()
        else
                n=sqrt(size(wfall_p_0,1));

                chern_p=zeros(1,level);
                chern_m=zeros(1,level);

                l_index=1:Nai
                wfall_p=wfall_p_0(:,l_index,:);
                wfall_m=wfall_m_0(:,l_index,:);
                umap_p=reshape(wfall_p,[n,n,size(wfall_p,3)]);
                umap_m=reshape(wfall_m,[n,n,size(wfall_p,3)]);

                nindex=1:n;
                nindex_inc=mod(nindex,n)+1;
                u1_p=umap_p(nindex,nindex,:);
                u2_p=umap_p(nindex_inc,nindex,:);
                u3_p=umap_p(nindex_inc,nindex_inc,:);
                u4_p=umap_p(nindex,nindex_inc,:);

                bcmap_p=angle(dot(u1_p,u2_p,3)...
                        .*dot(u2_p,u3_p,3)...
                        .*dot(u3_p,u4_p,3)...
                        .*dot(u4_p,u1_p,3));

                u1_m=umap_m(nindex,nindex,:);
                u2_m=umap_m(nindex_inc,nindex,:);
                u3_m=umap_m(nindex_inc,nindex_inc,:);
                u4_m=umap_m(nindex,nindex_inc,:);
                bcmap_m=angle(dot(u1_m,u2_m,3)...
                        .*dot(u2_m,u3_m,3)...
                        .*dot(u3_m,u4_m,3)...
                        .*dot(u4_m,u1_m,3));

                chern_p(l_index)=sum(bcmap_p(:))/(2*pi);
                chern_m(l_index)=sum(bcmap_m(:))/(2*pi);
        end
end