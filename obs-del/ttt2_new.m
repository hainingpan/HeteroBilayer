function z_final=ttt2(x,y,i,j,k,l)
%tensor product preserving certain index k,l, sum over i,j
%output index: first # of k follows the order of x, then goes the unsummed
%indices from x to y
if ~(isa(x,'tensor') || isa(x,'sptensor'))
    x=tensor(x);
end
if ~(isa(y,'tensor') || isa(y,'sptensor'))
    y=tensor(y);
end

if ~isempty(i) && ~isempty(j)
    assert(isequal(x.size(i),y.size(j)),sprintf('Must tensor product along the same dimensions #1 (%d)!= #2 (%d)',x.size(i),y.size(j)))
else
    if ~(isempty(i) && isempty(j))
        error('Must tensor product along the same dimensions');
    end
    assert(isempty(i) && isempty(j),sprintf('Must both be empty'))
end
assert(isempty(intersect(i, k)),sprintf('dim of #1 has overlap (%d,$d)',i,k));
assert(isempty(intersect(j, l)),sprintf('dim of #2 has overlap (%d,$d)',j,l));
assert(isequal(x.size(k),y.size(l)),sprintf('Must perserve along the same dimensions'));

nonpreserving_x = setdiff(1:ndims(x),k);
remdims_x=setdiff(nonpreserving_x,i);
nonpreserving_y = setdiff(1:ndims(y),l);
remdims_y=setdiff(nonpreserving_y,j);

xsize=x.size;
ysize=y.size;
preserve_size=xsize(k);

x = permute(x, [k,i,remdims_x]);
y = permute(y, [l,j,remdims_y]);

x_flatten=reshape(x,[prod(preserve_size),xsize([i,remdims_x])]);
y_flatten=reshape(y,[prod(preserve_size),ysize([j,remdims_y])]);
z_flatten=tenzeros([prod(preserve_size),remdims_x,remdims_y]);
for preserve_size_index=1:prod(preserve_size)
    x_slice_index=repmat({':'},1,ndims(x_flatten));
    x_slice_index{1}=preserve_size_index;
    x_slice=x_flatten(x_slice_index{:});
    y_slice_index=repmat({':'},1,ndims(y_flatten));
    y_slice_index{1}=preserve_size_index;
    y_slice=y_flatten(y_slice_index{:});
    z_slice=ttt(x_slice,y_slice,1:length(i),1:length(j)); % xsize(remdims_x), ysize(remdims_y)
    z_flatten_index=repmat({':'},1,ndims(z_flatten));
    z_flatten_index{1}=preserve_size_index;
    z_flatten(z_flatten_index{:})=z_slice;
end
z_final=reshape(z_flatten,preserve_size);
end