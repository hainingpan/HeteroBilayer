function [chern_p,chern_m]=chern_gs(wfall,valley_index)
% ignore q for now. ignore nonabelian wilson loop

n=sqrt(size(wfall,1));
NL=size(wfall,2);

% valley_p=repmat(valley_index==1,[1,1,size(wfall,3)]);
% valley_m=repmat(valley_index==-1,[1,1,size(wfall,3)]);

wfall_p=zeros(size(wfall,1),size(wfall,3)/2);
wfall_m=zeros(size(wfall,1),size(wfall,3)/2);
for k_index=1:size(wfall,1)
    for l_index=1:2
        if valley_index(k_index,l_index)==1
            wfall_p(k_index,:)=squeeze(wfall(k_index,l_index,1:end/2));
        elseif valley_index(k_index,l_index)==-1
            wfall_m(k_index,:)=squeeze(wfall(k_index,l_index,end/2+1:end));
        end
    end
end

umap_p=reshape(wfall_p,[n,n,size(wfall,3)/2]);
umap_m=reshape(wfall_m,[n,n,size(wfall,3)/2]);

nindex=1:n-1;
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

% bcmap_p=angle(dot(umap_p(1:end-1,1:end-1,:),umap_p(1:end-1,2:end,:),3)...
% .*dot(umap_p(1:end-1,2:end,:),umap_p(2:end,2:end,:),3)...
% .*dot(umap_p(2:end,2:end,:),umap_p(2:end,1:end-1,:),3)...
% .*dot(umap_p(2:end,1:end-1,:),umap_p(1:end-1,1:end-1,:),3));
% 
% bcmap_m=angle(dot(umap_m(1:end-1,1:end-1,:),umap_m(1:end-1,2:end,:),3)...
% .*dot(umap_m(1:end-1,2:end,:),umap_m(2:end,2:end,:),3)...
% .*dot(umap_m(2:end,2:end,:),umap_m(2:end,1:end-1,:),3)...
% .*dot(umap_m(2:end,1:end-1,:),umap_m(1:end-1,1:end-1,:),3));

chern_p=sum(bcmap_p(:))/(2*pi);
chern_m=sum(bcmap_m(:))/(2*pi);
end