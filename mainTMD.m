function params=mainTMD(varargin)
    p=inputParser;
    addParameter(p,'a_b',3.575e-10*5.076e6);
    addParameter(p,'a_t',3.32e-10*5.076e6);
    addParameter(p,'theta',0);
    addParameter(p,'m_b',0.65);
    addParameter(p,'m_t',0.35);
    addParameter(p,'V_b',4.1);
    addParameter(p,'V_t',0);
    addParameter(p,'psi_b',14);
    addParameter(p,'psi_t',0);
    addParameter(p,'w',1.3);
    addParameter(p,'Vz_b',0);
    addParameter(p,'Vz_t',-36.8);
    addParameter(p,'Nmax',5);
    addParameter(p,'omega',exp(1i*2*pi/3));
    addParameter(p,'valley',+1);% -1: -K; +1; +k
    addParameter(p,'nu',0);%filling factor number(#1) per site(#2); no spin degeneracy considered
    addParameter(p,'hole',1); %1: hole-like energy band ; -1: particle-like energy band
    addParameter(p,'n',15); % n*n grid in momentum space
    addParameter(p,'d',5); % gate to sample distance (nm)
    addParameter(p,'epsilon',1); % gate to sample distance
    addParameter(p,'tsymm',0); % T-symm is enforced 
    addParameter(p,'shift',0); % shift to Kb
    addParameter(p,'SDW',20e-3); % SDW
    
    parse(p,varargin{:});
    params=struct('a_b',p.Results.a_b,'a_t',p.Results.a_t,'theta',p.Results.theta*pi/180,'m_b',p.Results.m_b*0.511e6,'m_t',p.Results.m_t*0.511e6,'V_b',p.Results.V_b*1e-3,'V_t',p.Results.V_t*1e-3,'psi_b',p.Results.psi_b/360*2*pi,'psi_t',p.Results.psi_t/360*2*pi,'w',p.Results.w*1e-3,'Vz_b',p.Results.Vz_b*1e-3,'Vz_t',p.Results.Vz_t*1e-3,'Nmax',p.Results.Nmax,'omega',p.Results.omega,'valley',p.Results.valley,'nu',p.Results.nu,'hole',p.Results.hole,'n',p.Results.n,'d',p.Results.d*1e-9*5.076e6,'epsilon',p.Results.epsilon,'tsymm',p.Results.tsymm,'shift',p.Results.shift,'SDW',p.Results.SDW);

    delta=(params.a_b-params.a_t)/params.a_t;
    params.aM=params.a_b/sqrt(delta^2+params.theta^2);
    params.aM1=params.aM*[1,0];
    params.aM2=params.aM*[1/2,sqrt(3)/2];
    params.g1=[0,4*pi/(sqrt(3)*params.aM)];
    params.g2=[-2*pi/params.aM,(2*pi)/(sqrt(3)*params.aM)];
    params.g3=[-2*pi/params.aM,-(2*pi)/(sqrt(3)*params.aM)];
    params.g4=[0,-4*pi/(sqrt(3)*params.aM)];
    params.g5=[2*pi/params.aM,-(2*pi)/(sqrt(3)*params.aM)];
    params.g6=[2*pi/params.aM,(2*pi)/(sqrt(3)*params.aM)];
    params.bM1=params.g5;
    params.bM2=params.g1;
    neighbor_index=generate_shell(params.Nmax);
    params.area=sqrt(3)/2*params.aM^2;

    params.Kb=4*pi/(3*params.aM)*[-1/2,sqrt(3)/2];
    params.Kt=4*pi/(3*params.aM)*[1/2,sqrt(3)/2];

    if mod(params.n,3)~=0
        warning('n={%d} is not multiple of 3, which does not go through K point',params.n);
    end

    rotate=@(x) [cos(x) -sin(x);sin(x) cos(x)]; %rotate anticlockwise

    params.valley_polarized=0;    % unpolarized filling initially
    params.fermisurface=1;  % filled by Fermi surface initially
    params.auto_generate_q=1;
    params.span='b';
    params.chern='c';
    params.Sq_index=[];
    params.nq_p_index=[];
    params.nq_m_index=[];
    params.sxy=1; % magnitude of spin of xy component
    params.sz_p=0; %Canted AF
    params.sz_m=0; %Canted AF
    params.s0=0;
    params.phi_p=0; % sum cos(q.r+phi) for sigma_0 at +K
    params.phi_m=0; % sum cos(q.r+phi) for sigma_0 at -K
    params.NL=0;

    % for single-particle
    if params.nu==0  
        [ux,uy]=ndgrid(1:params.n,1:params.n);
        params.K_index=[(ux(:)-1)/(params.n),(uy(:)-1)/(params.n)];
        params.K=params.K_index*[params.bM1;params.bM2];
        ailist=[0,0];
        am_index=eye(2);
    end


    % QAHE(FM_z), spanned by b
    if params.nu==[1,1]
        ailist=[0,0];
        am_index=eye(2);
        params.valley_polarized=1;
        params.fermisurface=0;
        
        params.SDW=0;
        % params.sz=1;
    end
    % QAHE(FM_z) WC, cascaded
    if params.nu==[2,2]
        ailist=[[0,0];[1,0];[2,0]];
        am_index=[[1,1];[2,-1]]; % am=am_index* [aM1;aM2]; am_index=[am1_index,am2_index];
        params.valley_polarized=1;
        params.fermisurface=0;
        params.SDW=0;
    end
    % QAHE(FM_z) WC, spanned by q
    if params.nu==[4,4]
        ailist=[[0,0];[1,0];[2,0]];
        am_index=[[1,1];[2,-1]]; % am=am_index* [aM1;aM2]; am_index=[am1_index,am2_index];
        params.valley_polarized=1;
        params.fermisurface=0;
        params.span='q';
        params.auto_generate_q=0;
        params.SDW=0;
        % params.sz=1;
    end

    % FM_x without Wigner Crystal, spanned by b
    if params.nu==[-1,-1]
        ailist=[0,0];
        am_index=eye(2);
        % params.SDW=10e-3;
        q0_index=[0,0];
        params.Sq_index=[q0_index];
    end
    
    % FM_x  Wigner Crystal, cascaded
    if params.nu==[-2,-2] 
        ailist=[[0,0];[1,0];[2,0]];
        am_index=[[1,1];[2,-1]]; % am=am_index* [aM1;aM2]; am_index=[am1_index,am2_index];
        % params.SDW=10e-3;   % the strength of SDW 
        q0_index=[0,0];
        params.Sq_index=[q0_index];
    end

    % FM_x  Wigner Crystal, spanned by q
    if params.nu==[-4,-4] 
        ailist=[[0,0];[1,0];[2,0]];
        am_index=[[1,1];[2,-1]]; % am=am_index* [aM1;aM2]; am_index=[am1_index,am2_index];
        % params.SDW=10e-3;   % the strength of SDW 
        params.span='q';
        q0_index=[0,0];
        params.Sq_index=[q0_index];
        params.auto_generate_q=0;
    end

    %trivial Mott insulator, 120 AF, +1 chirality
    if params.nu==[3,3] 
        ailist=[[0,0];[1,0];[2,0]];
        am_index=[[1,1];[2,-1]]; % am=am_index* [aM1;aM2]; am_index=[am1_index,am2_index];
        % params.SDW=10e-3;   % the strength of SDW 
        Q=4*pi/(3*params.aM)*[1,0];
        q0=[Q;Q*rotate(2*pi/3);Q*rotate(4*pi/3)];
        params.Sq_index=q0/[params.bM1;params.bM2]; % SDW S(r)*tau, S(r)=sum_{q} {cos(q*r);sin(q*r)}
    end

    %trivial Mott insulator, 120 AF, spanned by q, +1 chirality
    if params.nu==[6,6] 
        ailist=[[0,0];[1,0];[2,0]];
        am_index=[[1,1];[2,-1]]; % am=am_index* [aM1;aM2]; am_index=[am1_index,am2_index];
        % params.SDW=10e-3;   % the strength of SDW 
        Q=4*pi/(3*params.aM)*[1,0];
        q0=[Q;Q*rotate(2*pi/3);Q*rotate(4*pi/3)];
        params.Sq_index=q0/[params.bM1;params.bM2]; % SDW S(r)*tau, S(r)=sum_{q} {cos(q*r);sin(q*r)}
        params.sz=0e-3;
        params.span='q';
        params.auto_generate_q=0;
    end

    %trivial Mott insulator, 120 AF, -1 chirality
    if params.nu==[-3,-3] 
        ailist=[[0,0];[1,0];[2,0]];
        am_index=[[1,1];[2,-1]]; % am=am_index* [aM1;aM2]; am_index=[am1_index,am2_index];
        % params.SDW=10e-3;   % the strength of SDW 
        Q=4*pi/(3*params.aM)*[-1,0];
        q0=[Q;Q*rotate(2*pi/3);Q*rotate(4*pi/3)];
        params.Sq_index=q0/[params.bM1;params.bM2]; % SDW S(r)*tau, S(r)=sum_{q} {cos(q*r);sin(q*r)}
    end
    %trivial Mott insulator, 120 AF, spanned by q, -1 chirality
    if params.nu==[-6,-6] 
        ailist=[[0,0];[1,0];[2,0]];
        am_index=[[1,1];[2,-1]]; % am=am_index* [aM1;aM2]; am_index=[am1_index,am2_index];
        % params.SDW=10e-3;   % the strength of SDW 
        Q=4*pi/(3*params.aM)*[-1,0];
        q0=[Q;Q*rotate(2*pi/3);Q*rotate(4*pi/3)];
        params.sz=0e-3;
        params.Sq_index=q0/[params.bM1;params.bM2]; % SDW S(r)*tau, S(r)=sum_{q} {cos(q*r);sin(q*r)}
        params.span='q';
        params.auto_generate_q=0;
    end

    % QSHE
    if params.nu==[2,1]
        ailist=[0,0];
        params.NL=1;
        am_index=eye(2);
        params.tsymm=1;
        params.SDW=0;
    end

    % Chern insulator
    if params.nu==[4,2]
        ailist=[0,0];
        params.NL=1;
        am_index=eye(2);
        % params.tsymm=1;
        params.nq_p_index=[[0,0]];
        params.nq_m_index=[[0,0]];
        params.s0=40e-3/params.SDW;
        params.sz_p=40e-3/params.SDW;
        params.sz_m=40e-3/params.SDW;
    end

    % Honeycomb FM
    if params.nu==[2,3]
        ailist=[[0,0];[1,0];[2,0]];
        am_index=[[1,1];[2,-1]]; % am=am_index* [aM1;aM2]; am_index=[am1_index,am2_index];
        Q=4*pi/(3*params.aM)*[1,0];
        q0=[Q;Q*rotate(2*pi/3);Q*rotate(4*pi/3)];
        params.nq_p_index=q0/[params.bM1;params.bM2]; 
        params.phi_p=pi;
        params.nq_m_index=q0/[params.bM1;params.bM2]; 
        params.phi_m=pi;
        params.span='q';
        params.auto_generate_q=0;
        params.SDW=10e-3;
        params.Sq_index=[];
        params.s0=15e-3/params.SDW;
        params.sz_p=1;
        params.sz_m=1;
    end
    
    % Honeycomb AFM
    if params.nu==[-2,-3]
        ailist=[[0,0];[1,0];[2,0]];
        am_index=[[1,1];[2,-1]]; % am=am_index* [aM1;aM2]; am_index=[am1_index,am2_index];
        Q=4*pi/(3*params.aM)*[1,0];
        q0=[Q;Q*rotate(2*pi/3);Q*rotate(4*pi/3)];
        params.nq_p_index=q0/[params.bM1;params.bM2]; 
        params.phi_p=-2*pi/3;
        params.nq_m_index=q0/[params.bM1;params.bM2]; 
        params.phi_m=2*pi/3;

        params.span='q';
        params.auto_generate_q=0;
        params.SDW=10e-3;
        params.Sq_index=[];
        params.s0=10e-3/params.SDW;
        params.sz_p=0;
        params.sz_m=0;
        params.NL=1;
    end

    % Metalic state 
    if params.nu==[4,6]
        ailist=[[0,0];[1,0];[2,0]];
        am_index=[[1,1];[2,-1]]; % am=am_index* [aM1;aM2]; am_index=[am1_index,am2_index];
        params.valley_polarized=1;
        params.fermisurface=0;
        params.span='q';
        params.auto_generate_q=0;
        params.SDW=0;
        % params.sz=1;
    end

    % h1 bM1+ h2 bM2
    h1index=neighbor_index(:,1);
    h2index=neighbor_index(:,2);
    params.h1index=h1index;
    params.h2index=h2index;
    [h1matX,h1matY]=ndgrid(h1index(:));
    [h2matX,h2matY]=ndgrid(h2index(:));

    params.DeltaTmat=DeltaT(h1matX-h1matY,h2matX-h2matY,params);
    params.DeltaTTmat=DeltaTT(h1matX-h1matY,h2matX-h2matY,params);
    params.Deltabmat=Deltal(h1matX-h1matY,h2matX-h2matY,1,params);
    params.Deltatmat=Deltal(h1matX-h1matY,h2matX-h2matY,-1,params);
    
    valley0=params.valley;
    R=[[0,-1];[1,0]];
    params.ailist=ailist;
    bm_index=1/size(ailist,1)*R'*am_index*R;  %check oneNotes
    params.am1=am_index(1,:)*[params.aM1;params.aM2];
    params.am2=am_index(2,:)*[params.aM1;params.aM2];
    params.bm1=bm_index(1,:)*[params.bM1;params.bM2];
    params.bm2=bm_index(2,:)*[params.bM1;params.bM2];
    
    params.kb=-1/3*params.bm1+1/3*params.bm2;
    params.kt=1/3*params.bm1+2/3*params.bm2;
    m=(params.kb+params.kt)/2;
    gamma=[0,0];
    kt_m_x=linspace(params.kt(1),m(1),40);
    kt_m_y=linspace(params.kt(2),m(2),40);
    m_kb_x=linspace(m(1),params.kb(1),40);
    m_kb_y=linspace(m(2),params.kb(2),40);
    kb_gamma_x=linspace(params.kb(1),gamma(1),80);
    kb_gamma_y=linspace(params.kb(2),gamma(2),80);
    k_line_x=[kt_m_x,m_kb_x(2:end),kb_gamma_x(2:end)];
    k_line_y=[kt_m_y,m_kb_y(2:end),kb_gamma_y(2:end)];
    params.k_line=[k_line_x(:),k_line_y(:)];

    kt0=params.kb;
    kt1=kt0*rotate(pi*2/3);
    kt2=kt0*rotate(-pi*2/3);
    kt0_kt1_x=linspace(kt0(1),kt1(1),40);
    kt0_kt1_y=linspace(kt0(2),kt1(2),40);
    kt1_kt2_x=linspace(kt1(1),kt2(1),40);
    kt1_kt2_y=linspace(kt1(2),kt2(2),40);
    kt2_kt0_x=linspace(kt2(1),kt0(1),40);
    kt2_kt0_y=linspace(kt2(2),kt0(2),40);
    k_line_C3_x=[kt0_kt1_x,kt1_kt2_x(2:end),kt2_kt0_x(2:end)];
    k_line_C3_y=[kt0_kt1_y,kt1_kt2_y(2:end),kt2_kt0_y(2:end)];
    params.k_line_C3=[k_line_C3_x(:),k_line_C3_y(:)];

    
    if params.nu~=0
        if params.NL==0
            params.NL=size(params.ailist,1)*params.nu(1)/params.nu(2);
        end

        if params.auto_generate_q==1
            %square [0,1]x[0,1] transformed to new shape
            new_pts=int8([[0,0];[1,0];[1,1];[0,1]]/bm_index);
            %check whether inside the transformed square
            xrange=min(new_pts(:,1)):max(new_pts(:,1));
            yrange=min(new_pts(:,2)):max(new_pts(:,2));
            [xgrid,ygrid]=ndgrid(xrange,yrange);
            [in_sq,on_sq]=inpolygon(xgrid(:), ygrid(:), new_pts(:,1), new_pts(:,2));
            in_sq=in_sq&(~on_sq);
            q_index_sq=[xgrid(in_sq),ygrid(in_sq)];
            %check whether on the line [0,0]->[1,0]
            if new_pts(2,1)-1>0
                xlist=int8(linspace(1,new_pts(2,1)-1,new_pts(2,1)-1));
                xlist_rem=mod(xlist*new_pts(2,2),new_pts(2,1));
                in_line1_x=xlist_rem(xlist_rem==0);
                in_line1_y=in_line1_x*new_pts(2,2)/new_pts(2,1);
                q_index_line1=[in_line1_x(:),in_line1_y(:)];
            else
                q_index_line1=[];
            end
            %check whether on the line [0,0]->[0,1]
            if new_pts(3,1)-1>0
                xlist=int8(linspace(1,new_pts(3,1)-1,new_pts(3,1)-1));
                xlist_rem=mod(xlist*new_pts(3,2),new_pts(3,1));
                in_line2_x=xlist_rem(xlist_rem==0);
                in_line2_y=in_line2_x*new_pts(3,2)/new_pts(3,1);
                q_index_line2=[in_line2_x(:),in_line2_y(:)];
            else
                q_index_line2=[];
            end
            q_index=[0,0;q_index_sq;q_index_line1;q_index_line2];
            assert(size(q_index,1)==size(ailist,1),'q_index and ailist are not equal length');
            q_index=double(q_index)*bm_index;
        
            params.q_index=q_index;
        else
            params.q_index=neighbor_index*bm_index;
        end
        params.q=params.q_index*[params.bM1;params.bM2]; %% !!! This may be a problem when bm1 and bm2 are not separeted by 120 deg; in such case, it needs to generate a neighbor_index by inputing the two reciprocal lattice vectors.
        Nq=size(params.q,1);

        if params.span=='b'
            params.b_index=neighbor_index;
        elseif params.span=='q'
            params.b_index=[[0,0]];
        end
        params.b=params.b_index*[params.bM1;params.bM2];
        Nb=size(params.b,1);

        [ux,uy]=ndgrid(1:params.n,1:params.n);
        params.k_index=[(2*ux(:)-params.n-1)/(2*params.n),(2*uy(:)-params.n-1)/(2*params.n)];
        % params.k_index=[(2*ux(:)-params.n)/(2*params.n),(2*uy(:)-params.n)/(2*params.n)];
        % params.k_index=[(ux(:)-1)/(params.n),(uy(:)-1)/(params.n)];

        params.k=params.k_index*[params.bm1;params.bm2];

        [ux,uy]=ndgrid(1:2*params.n,1:2*params.n);
        params.k_dense_index=[(ux(:)-1)/(2*params.n),(uy(:)-1)/(2*params.n)];
        params.k_dense=params.k_dense_index*[params.bm1;params.bm2];

        [ux,uy]=ndgrid(0:params.n,0:params.n);
        params.k_index_bc=[(2*ux(:)-params.n)/(2*params.n),(2*uy(:)-params.n)/(2*params.n)];
        % params.k_index_bc=[(ux(:)-1)/(params.n),(uy(:)-1)/(params.n)];
        params.k_bc=params.k_index_bc*[params.bm1;params.bm2];

        [q_a_x,b_a_x,q_b_x,b_b_x]=ndgrid(params.q_index(:,1),params.b_index(:,1),params.q_index(:,1),params.b_index(:,1));
        [q_a_y,b_a_y,q_b_y,b_b_y]=ndgrid(params.q_index(:,2),params.b_index(:,2),params.q_index(:,2),params.b_index(:,2));
        params.valley=1;
        h1=reshape(b_a_x+q_a_x-b_b_x-q_b_x,[Nq*Nb,Nq*Nb]);
        h2=reshape(b_a_y+q_a_y-b_b_y-q_b_y,[Nq*Nb,Nq*Nb]);
        params.Delta_T_p=DeltaT(h1,h2,params);
        params.Delta_TT_p=DeltaTT(h1,h2,params);
        params.Delta_b_p=Deltal(h1,h2,1,params);
        params.Delta_t_p=Deltal(h1,h2,-1,params);
        params.valley=-1;
        params.Delta_T_m=DeltaT(h1,h2,params);
        params.Delta_TT_m=DeltaTT(h1,h2,params);
        params.Delta_b_m=Deltal(h1,h2,1,params);
        params.Delta_t_m=Deltal(h1,h2,-1,params);
        params.valley=valley0;

        if params.SDW~=0
            [q_a_x,b_a_x,q_b_x,b_b_x]=ndgrid(params.q_index(:,1),params.b_index(:,1),params.q_index(:,1),params.b_index(:,1));
            [q_a_y,b_a_y,q_b_y,b_b_y]=ndgrid(params.q_index(:,2),params.b_index(:,2),params.q_index(:,2),params.b_index(:,2));
            h1=b_a_x+q_a_x-b_b_x-q_b_x;
            h2=b_a_y+q_a_y-b_b_y-q_b_y;
            S_pm=kron([1,0;0,0],reshape(S(h1(:),h2(:),1,-1,params),[Nq*Nb,Nq*Nb]));
            S_mp=kron([1,0;0,0],reshape(S(h1(:),h2(:),-1,1,params),[Nq*Nb,Nq*Nb]));
            S_pp=kron([1,0;0,0],reshape(S(h1(:),h2(:),1,1,params),[Nq*Nb,Nq*Nb]));
            S_mm=kron([1,0;0,0],reshape(S(h1(:),h2(:),-1,-1,params),[Nq*Nb,Nq*Nb]));
            params.S_tau=[S_pp,S_pm;S_mp,S_mm];
        else
            params.S_tau=0;
        end

        % For the detailed def, check oneNotes
        [q_a_x,q_b_x,q_g_x,q_d_x,b_a_x,b_b_x,b_g_x,b_d_x]=ndgrid(params.q_index(:,1),params.q_index(:,1),params.q_index(:,1),params.q_index(:,1),params.b_index(:,1),params.b_index(:,1),params.b_index(:,1),params.b_index(:,1));
        delta_x=-q_a_x-b_a_x+q_b_x+b_b_x+q_g_x+b_g_x-q_d_x-b_d_x;
        [q_a_y,q_b_y,q_g_y,q_d_y,b_a_y,b_b_y,b_g_y,b_d_y]=ndgrid(params.q_index(:,2),params.q_index(:,2),params.q_index(:,2),params.q_index(:,2),params.b_index(:,2),params.b_index(:,2),params.b_index(:,2),params.b_index(:,2));
        delta_y=-q_a_y-b_a_y+q_b_y+b_b_y+q_g_y+b_g_y-q_d_y-b_d_y;
        delta_tensor1=int8((abs(delta_x)<1e-10) & (abs(delta_y)<1e-10)); %-q_a-b_a+q_b+b_b+q_g+b_g-q_d-b_d, index: q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d
        clear delta_x delta_y
        params.delta_tensor1=tensor(delta_tensor1,[Nq,Nq,Nq,Nq,Nb,Nb,Nb,Nb]);
        delta_x=-q_a_x-b_a_x+q_b_x+b_b_x-q_g_x-b_g_x+q_d_x+b_d_x;
        clear q_a_x b_a_x q_b_x b_b_x q_g_x b_g_x q_d_x b_d_x
        delta_y=-q_a_y-b_a_y+q_b_y+b_b_y-q_g_y-b_g_y+q_d_y+b_d_y;
        clear q_a_y b_a_y q_b_y b_b_y q_g_y b_g_y q_d_y b_d_y
        delta_tensor2=int8((abs(delta_x)<1e-10) & (abs(delta_y)<1e-10)); %-q_a-b_a+q_b+b_b-q_g-b_g+q_d+b_d, index: q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d
        params.delta_tensor2=tensor(delta_tensor2,[Nq,Nq,Nq,Nq,Nb,Nb,Nb,Nb]);
        alpha=0.00729735; % eV*nm

        [q_g_x,q_d_x,b_g_x,b_d_x]=ndgrid(params.q(:,1),params.q(:,1),params.b(:,1),params.b(:,1));
        [q_g_y,q_d_y,b_g_y,b_d_y]=ndgrid(params.q(:,2),params.q(:,2),params.b(:,2),params.b(:,2));
        qd=params.d*sqrt((b_g_x+q_g_x-b_d_x-q_d_x).^2+(b_g_y+q_g_y-b_d_y-q_d_y).^2);
        params.V1=alpha*2*pi*tanh(qd+1e-18)./(qd+1e-18)*params.d;
        assert(sum(isnan(params.V1(:)))==0,'V1 has NaN');
        params.V1=tensor(params.V1,[Nq,Nq,Nb,Nb]);
        Nk=size(params.k,1);
        [k_a_x,k_b_x,q_a_x,q_d_x,b_a_x,b_d_x]=ndgrid(params.k(:,1),params.k(:,1),params.q(:,1),params.q(:,1),params.b(:,1),params.b(:,1));
        k_x=k_a_x+b_d_x+q_d_x-k_b_x-b_a_x-q_a_x;
        clear k_a_x b_d_x q_d_x k_b_x b_a_x q_a_x 
        [k_a_y,k_b_y,q_a_y,q_d_y,b_a_y,b_d_y]=ndgrid(params.k(:,2),params.k(:,2),params.q(:,2),params.q(:,2),params.b(:,2),params.b(:,2));
        k_y=k_a_y+b_d_y+q_d_y-k_b_y-b_a_y-q_a_y;
        clear k_a_y b_d_y q_d_y k_b_y b_a_y q_a_y
        qd=params.d*sqrt(k_x.^2+k_y.^2);
        clear k_x k_y
        params.V2=alpha*2*pi*tanh(qd+1e-18)./(qd+1e-18)*params.d;
        clear qd;
        assert(sum(isnan(params.V2(:)))==0,'V2 has NaN');
        params.V2=tensor(params.V2,[Nk,Nk,Nq,Nq,Nb,Nb]);
    end