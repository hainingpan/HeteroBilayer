function re=S(h1,h2,tau1,tau2,params)
    % tau==+1 : +K
    % tau==-1 : -K

    sum=0*h1;

    if tau1*tau2==-1
        for i=1:size(params.Sq_index,1)
            sum=sum+iseq(h1,-params.Sq_index(i,1)*tau1).*iseq(h2,-params.Sq_index(i,2)*tau1);
        end
        sum=params.sxy*sum;
    else
        % for i=1:size(params.nq_index,1)
        %     sum=sum+iseq(h1,params.nq_index(i,1)).*iseq(h2,params.nq_index(i,2))*exp(1i*params.phi)+iseq(h1,-params.nq_index(i,1)).*iseq(h2,-params.nq_index(i,2))*exp(-1i*params.phi);
        % end
        % sum=params.s0*sum/2+tau1*params.sz*(iseq(h1,0).*iseq(h2,0));
        if tau1==1
            for i=1:size(params.nq_p_index,1)
                sum=sum+iseq(h1,params.nq_p_index(i,1)).*iseq(h2,params.nq_p_index(i,2))*exp(1i*params.phi_p)+iseq(h1,-params.nq_p_index(i,1)).*iseq(h2,-params.nq_p_index(i,2))*exp(-1i*params.phi_p);
            end
            sum=params.s0*sum/2+params.sz_p*(iseq(h1,0).*iseq(h2,0));
        else
            for i=1:size(params.nq_m_index,1)
                sum=sum+iseq(h1,params.nq_m_index(i,1)).*iseq(h2,params.nq_m_index(i,2))*exp(1i*params.phi_m)+iseq(h1,-params.nq_m_index(i,1)).*iseq(h2,-params.nq_m_index(i,2))*exp(-1i*params.phi_m);
            end
            sum=params.s0*sum/2-params.sz_m*(iseq(h1,0).*iseq(h2,0));
        end
    end
    re=params.SDW*sum;
end