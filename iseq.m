function re=iseq(a,b)
% test whether a and b are equal within epsilon
% a: vector, b: scalar
    epsilon=1e-5;
    re=abs(a-b)<epsilon;
end

