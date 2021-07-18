function z_final=ttt2(x,y,i,j,k,l)
%tensor product preserving certain index k,l, sum over i,j
%output index: first # of k follows the order of x, then goes the unsummed
%indices from x to y
if isa(x,'double') || isa(x,'int8')
    x=tensor(x);
end

if isa(y,'double') || isa(y,'int8')
    y=tensor(y);
end

if ~isempty(i) & ~isempty(j)
    if x.size(i)~=y.size(j)
        error('Must tensor product along the same dimensions');
    end
else
    if ~(isempty(i) & isempty(j))
        error('Must tensor product along the same dimensions');
    end
end

if x.size(k)~=y.size(l)
    error('Must perserve along the same dimensions');
end
remdims_x = setdiff(1:ndims(x),k);
remdims_y = setdiff(1:ndims(y),l);
xsize=x.size;
ysize=y.size;

x = permute(x, [k,remdims_x]);
y = permute(y, [l,remdims_y]);

preserve_size=xsize(k);

x_mat=reshape(x.data,[prod(preserve_size),prod(xsize(remdims_x))]);
y_mat=reshape(y.data,[prod(preserve_size),prod(ysize(remdims_y))]);

% ind_x=find(ismember(remdims_x,i));
% ind_y=find(ismember(remdims_y,j));
ind_x=arrayfun(@(xx) find(remdims_x==xx) ,i);
ind_y=arrayfun(@(xx) find(remdims_y==xx) ,j);
tenprodsize=[xsize(setdiff(remdims_x,i)),ysize(setdiff(remdims_y,j))];

z=zeros([prod(preserve_size),prod(tenprodsize)]);

%easy matrix multiplication


for i=1:prod(preserve_size)
    if length(remdims_x)==1
        x_reshape=x_mat(i,:);
    else
        x_reshape=reshape(x_mat(i,:),xsize(remdims_x));
    end
    
    if length(remdims_y)==1
        y_reshape=y_mat(i,:);
    else
        y_reshape=reshape(y_mat(i,:),ysize(remdims_y));
    end
    if length(remdims_x)==1
        z_tensor=ttv(tensor(y_reshape),x_reshape',ind_y);
    end
    if length(remdims_y)==1
        z_tensor=ttv(tensor(x_reshape),y_reshape',ind_x);
    end
    if length(remdims_x)~=1 & length(remdims_y)~=1
        z_tensor=ttt(tensor(x_reshape),tensor(y_reshape),ind_x,ind_y);
    end
    z(i,:)=z_tensor(:);
end
z_final=reshape(z,[preserve_size,tenprodsize]);
