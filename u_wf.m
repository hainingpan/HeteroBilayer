function [u_b_p_grid,u_t_p_grid,u_b_m_grid,u_t_m_grid]=u_wf(vec,rx,ry,parameters)
    %rx,ry: scalar|array real space position
    [~,Nstate]=size(vec);
    [Nx,Ny]=size(rx);
    h1index=parameters.h1index;
    h2index=parameters.h2index;

    expo=exp(1i*(h1index(:)*parameters.bM1+h2index(:)*parameters.bM2)*[rx(:),ry(:)]');
    vec_b_p=vec(1:size(parameters.h1index,1),:);
    vec_t_p=vec(size(parameters.h1index,1)+1:2*size(parameters.h1index,1),:);
    vec_b_m=vec(2*size(parameters.h1index,1)+1:3*size(parameters.h1index,1),:);
    vec_t_m=vec(3*size(parameters.h1index,1)+1:4*size(parameters.h1index,1),:);
    u_b_p=vec_b_p.'*expo;
    u_t_p=vec_t_p.'*expo;
    u_b_m=vec_b_m.'*expo;
    u_t_m=vec_t_m.'*expo;
    
    u_b_p_grid=reshape(u_b_p,Nstate,Nx,Ny);
    u_t_p_grid=reshape(u_t_p,Nstate,Nx,Ny);    %Normalization to V_UC
    u_b_m_grid=reshape(u_b_m,Nstate,Nx,Ny);
    u_t_m_grid=reshape(u_t_m,Nstate,Nx,Ny);    %Normalization to V_UC
    
    u_b_p_grid=permute(u_b_p_grid,[2,3,1]);
    u_t_p_grid=permute(u_t_p_grid,[2,3,1]);
    u_b_m_grid=permute(u_b_m_grid,[2,3,1]);
    u_t_m_grid=permute(u_t_m_grid,[2,3,1]);
    
    end