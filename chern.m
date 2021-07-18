function ch=chern(wfall,valley_index,parameters)
    n=sqrt(size(wfall,1));
    NL=size(wfall,2);
    % level=1:NL*parameters.nu(1)/(2*parameters.nu(2));
    level=1:2;
    valley_p=repmat(valley_index==1,[1,1,size(wfall,3)]);
    valley_m=repmat(valley_index==-1,[1,1,size(wfall,3)]);
    wfall_p=wfall(valley_p);
    wfall_m=wfall(valley_m);
    wfall_p=reshape(wfall_p,[n,n,NL/2,size(wfall,3)]);
    wfall_m=reshape(wfall_m,[n,n,NL/2,size(wfall,3)]);
    umap=cat(3,wfall_p(:,:,1,:),wfall_m(:,:,1,:));
    % umap=reshape(wfall_2,[n,n,2,size(wfall,3)]); %Nkx,Nky,level index, wf component
    % umap=(wfmap(:,:,level,:));
    
    berrycur=zeros(n-1);
    for i=1:n
        for j=1:n
            u1=reshape(umap(i,j,:,:),[level(end),NL]);        
    %         u2=squeeze(umap(mod(i,n)+1,j,:,:))*((i+1>n)*kron(eye(4),[0,1;1,0])+(i+1<=n)*eye(8));        
    %         u3=squeeze(umap(mod(i,n)+1,mod(j,n)+1,:,:))*((i+1>n)*kron(eye(4),[0,1;1,0])+(i+1<=n)*eye(8))*...
    %             ((j+1>n)*kron(eye(2),kron([0,1;1,0],eye(2)))+(j+1<=n)*eye(8));        
    %         u4=squeeze(umap(i,mod(j,n)+1,:,:))*((j+1>n)*kron(eye(2),kron([0,1;1,0],eye(2)))+(j+1<=n)*eye(8));
            u2=reshape(umap(mod(i,n)+1,j,:,:),[level(end),NL])*((i+1>n)*kron(eye(2),parameters.perm1)+(i+1<=n)*eye(NL));        
            u3=reshape(umap(mod(i,n)+1,mod(j,n)+1,:,:),[level(end),NL])*((i+1>n)*kron(eye(2),parameters.perm1)+(i+1<=n)*eye(NL))*...
                ((j+1>n)*kron(eye(2),parameters.perm2)+(j+1<=n)*eye(NL)); 
            u4=reshape(umap(i,mod(j,n)+1,:,:),[level(end),NL])*((j+1>n)*kron(eye(2),parameters.perm2)+(j+1<=n)*eye(NL));
    
            U12=conj(u1)*u2.';
            U23=conj(u2)*u3.';
            U34=conj(u3)*u4.';
            U41=conj(u4)*u1.';
            berrycur(i,j)=angle(det(U12*U23*U34*U41));
        end
    end
    
    ch=sum(berrycur(:))/(2*pi);
    end