function T0=T_gen(k_beta_set,shift,params)
b_set=params.b;
q_set=params.q;
Nb=size(b_set,1);
Nq=size(q_set,1);
Nk=size(k_beta_set,1);
kb=params.kb;
kt=params.kt;
m_b=params.m_b;
m_t=params.m_t;
T0=zeros(Nb*Nq*2*2,Nb*Nq*2*2,Nk);

for k_beta_index=1:Nk
    kx=k_beta_set(k_beta_index,1);
    ky=k_beta_set(k_beta_index,2);
    [qx_set_grid,bx_set_grid]=ndgrid(q_set(:,1),b_set(:,1));
    [qy_set_grid,by_set_grid]=ndgrid(q_set(:,2),b_set(:,2));
    if shift==1
        k_b_p_x=kx+qx_set_grid(:)+bx_set_grid(:);
        k_b_p_y=ky+qy_set_grid(:)+by_set_grid(:);

        k_b_m_x=kx+qx_set_grid(:)+bx_set_grid(:);
        k_b_m_y=ky+qy_set_grid(:)+by_set_grid(:);

        k_t_p_x=kx+qx_set_grid(:)+bx_set_grid(:)-(kt(1)-kb(1));
        k_t_p_y=ky+qy_set_grid(:)+by_set_grid(:)-(kt(2)-kb(2));

        k_t_m_x=kx+qx_set_grid(:)+bx_set_grid(:)+(kt(1)-kb(1));
        k_t_m_y=ky+qy_set_grid(:)+by_set_grid(:)+(kt(2)-kb(2));
    else
        k_b_p_x=kx+qx_set_grid(:)+bx_set_grid(:)-kb(1);
        k_b_p_y=ky+qy_set_grid(:)+by_set_grid(:)-kb(2);
    
        k_t_p_x=kx+qx_set_grid(:)+bx_set_grid(:)-kt(1);
        k_t_p_y=ky+qy_set_grid(:)+by_set_grid(:)-kt(2);
    
        k_b_m_x=-kx-qx_set_grid(:)-bx_set_grid(:)-kb(1);
        k_b_m_y=-ky-qy_set_grid(:)-by_set_grid(:)-kb(2);
    
        k_t_m_x=-kx-qx_set_grid(:)-bx_set_grid(:)-kt(1);
        k_t_m_y=-ky-qy_set_grid(:)-by_set_grid(:)-kt(2);
    end
    T_b_m=-1/(2*m_b)*diag(k_b_m_x.^2+k_b_m_y.^2);
    T_t_m=-1/(2*m_t)*diag(k_t_m_x.^2+k_t_m_y.^2);
    T_m=[T_b_m,0*T_b_m;0*T_b_m,T_t_m];
    T_b_p=-1/(2*m_b)*diag(k_b_p_x.^2+k_b_p_y.^2);
    T_t_p=-1/(2*m_t)*diag(k_t_p_x.^2+k_t_p_y.^2);
    T_p=[T_b_p,0*T_b_p;0*T_b_p,T_t_p];
    T_tau=[T_p,0*T_m;0*T_m,T_m];
    T0(:,:,k_beta_index)=T_tau;
end
