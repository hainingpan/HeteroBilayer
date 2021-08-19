function re=S(h1,h2,tau1,tau2,params)
    % tau==+1 : +K
    % tau==-1 : -K

    sum=0;
    if tau1*tau2==-1
        for i=1:size(params.Sq_index,1)
            sum=sum+iseq(h1,-params.Sq_index(i,1)*tau1).*iseq(h2,-params.Sq_index(i,2)*tau1);
        end
    else
        % sum=0*h1;
        sum=tau1*params.sz*eye(size(h1));
    end
    re=params.SDW*sum;
end