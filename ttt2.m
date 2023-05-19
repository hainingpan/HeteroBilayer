function z_final=ttt2(x,y,i,j,k,l)
%tensor product preserving certain index k,l, sum over i,j
%output index: first # of k follows the order of x, then goes the unsummed
%indices from x to y
if ~isa(x,'tensor')
    x=tensor(x);
end

if ~isa(y,'tensor')
    y=tensor(y);
end

if ~isempty(i) && ~isempty(j)
    assert(isequal(x.size(i),y.size(j)),sprintf('Must tensor product along the same dimensions #1 (%d)!= #2 (%d)',x.size(i),y.size(j)))
else
    assert(isempty(i) && isempty(j),sprintf('Must both be empty'))
end

assert(isempty(intersect(i, k)),sprintf('dim of #1 has overlap (%d,$d)',i,k));
assert(isempty(intersect(j, l)),sprintf('dim of #2 has overlap (%d,$d)',j,l));
assert(isequal(x.size(k),y.size(l)),sprintf('Must perserve along the same dimensions'));

remdims_x = setdiff(1:ndims(x),k);
remdims_y = setdiff(1:ndims(y),l);
xsize=x.size;
ysize=y.size;

x = permute(x, [k,remdims_x]);
y = permute(y, [l,remdims_y]);

preserve_size=xsize(k);

x_mat=reshape(x.data,[prod(preserve_size),prod(xsize(remdims_x))]);
y_mat=reshape(y.data,[prod(preserve_size),prod(ysize(remdims_y))]);

ind_x=arrayfun(@(xx) find(remdims_x==xx) ,i);
ind_y=arrayfun(@(xx) find(remdims_y==xx) ,j);
tenprodsize=[xsize(setdiff(remdims_x,i)),ysize(setdiff(remdims_y,j))];

z=zeros([prod(preserve_size),prod(tenprodsize)]);

%easy matrix multiplication


for i=1:prod(preserve_size)
    if length(remdims_x)==1
        x_reshape=x_mat(i,:);
    else
        x_reshape=tensor(x_mat(i,:),xsize(remdims_x));
    end
    
    if length(remdims_y)==1
        y_reshape=y_mat(i,:);
    else
        y_reshape=tensor(y_mat(i,:),ysize(remdims_y));
    end
    if length(remdims_x)==1
        z_tensor=ttv(tensor(y_reshape),x_reshape',ind_y);
    end
    if length(remdims_y)==1
        z_tensor=ttv(tensor(x_reshape),y_reshape',ind_x);
    end
    if length(remdims_x)~=1 & length(remdims_y)~=1
        z_tensor=ttt((x_reshape),(y_reshape),ind_x,ind_y);
    end
    z(i,:)=z_tensor(:);
end
% clear x_mat y_mat x_mat x_reshape y_reshape x y

if length([preserve_size,tenprodsize])==1 %the final result is simply a vector 
    z_final=z;
else
    z_final=tensor(z,[preserve_size,tenprodsize]);
end