function [s0,sx,sy,sz]=S_r(ave2,rmap_x,rmap_y,l,params)
    % l==1, bottom;
    % l==2, top;
    Nk=size(params.k,1);
    % Nq=size(params.q,1);
    Nai=size(params.ailist,1);  % The expansion of super cell
    [q_a_x,q_b_x,b_a_x,b_b_x,r_x]=ndgrid(params.q(:,1),params.q(:,1),params.b(:,1),params.b(:,1),rmap_x);
    [q_a_y,q_b_y,b_a_y,b_b_y,r_y]=ndgrid(params.q(:,2),params.q(:,2),params.b(:,2),params.b(:,2),rmap_y);
    expm=exp(-1i*((q_a_x-q_b_x+b_a_x-b_b_x).*r_x+(q_a_y-q_b_y+b_a_y-b_b_y).*r_y));    %q_a,q_b,b_a,b_b,r
    %ave2: k_a,q_d,b_d,l_a,t_a,q_g,b_g,l_b,t_b
    % size_ave2=size(ave2);
    % size_ave2([1,4,8])=[];
    ave2=collapse(ave2(:,:,:,l,:,:,:,l,:),[1]);  % q_d,b_d,t_a,q_g,b_g,t_b
    expm_ave2=ttt2(expm,ave2,[1,2,3,4],[1,4,2,5],[],[]); % r,t_a,t_b
    sigma0=eye(2);
    sigmax=[0,1;1,0];
    sigmay=[0,-1i;1i,0];
    sigmaz=[1,0;0,-1];
    
    s0=ttt2(sigma0,expm_ave2,[1,2],[2,3],[],[])/(Nk*Nai);
    sx=ttt2(sigmax,expm_ave2,[1,2],[2,3],[],[])/(Nk*Nai);
    sy=ttt2(sigmay,expm_ave2,[1,2],[2,3],[],[])/(Nk*Nai);
    sz=ttt2(sigmaz,expm_ave2,[1,2],[2,3],[],[])/(Nk*Nai);
    assert(max(abs(imag(s0(:))))<1e-10,sprintf('s0 not real %e',max(abs(imag(s0(:))))));
    s0=real(s0);
    assert(max(abs(imag(sx(:))))<1e-10,sprintf('sx not real %e',max(abs(imag(sx(:))))));
    sx=real(sx);
    assert(max(abs(imag(sy(:))))<1e-10,sprintf('sy not real %e',max(abs(imag(sy(:))))));
    sy=real(sy);
    assert(max(abs(imag(sz(:))))<1e-10,sprintf('sz not real %e',max(abs(imag(sz(:))))));
    sz=real(sz);
end