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
    addParameter(p,'valley',-1);% -1: -K; +1; +k
    addParameter(p,'nu',[1,1]);%filling factor number(#1) per site(#2); no spin degeneracy considered
    addParameter(p,'hole',1); %1: hole-like energy band ; -1: particle-like energy band
    
    parse(p,varargin{:});
    params=struct('a_b',p.Results.a_b,'a_t',p.Results.a_t,'theta',p.Results.theta*pi/180,'m_b',p.Results.m_b*0.511e6,'m_t',p.Results.m_t*0.511e6,'V_b',p.Results.V_b*1e-3,'V_t',p.Results.V_t*1e-3,'psi_b',p.Results.psi_b/360*2*pi,'psi_t',p.Results.psi_t/360*2*pi,'w',p.Results.w*1e-3,'Vz_b',p.Results.Vz_b*1e-3,'Vz_t',p.Results.Vz_t*1e-3,'Nmax',p.Results.Nmax,'omega',p.Results.omega,'valley',p.Results.valley,'nu',p.Results.nu,'hole',p.Results.hole);

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
    params.kb=4*pi/(3*params.aM)*[-1/2,sqrt(3)/2];
    params.kt=4*pi/(3*params.aM)*[1/2,sqrt(3)/2];
    % if params.valley==-1
    %     params.kb0=4*pi/(3*params.aM)*[1/2,sqrt(3)/2];
    % else
    %     params.kb0=4*pi/(3*params.aM)*[-1/2,sqrt(3)/2];
    % end
    % params.kb=params.kb-params.kb0;
    % params.kt=params.kt-params.kb0;

    % Nrange=-params.Nmax:params.Nmax;
%     [h1index,h2index]=meshgrid(Nrange,Nrange);
    neighbor_index=generate_shell(params.Nmax);
    h1index=neighbor_index(:,1);
    h2index=neighbor_index(:,2);
    params.h1index=h1index;
    params.h2index=h2index;
    [h1matX,h1matY]=meshgrid(h1index(:));
    [h2matX,h2matY]=meshgrid(h2index(:));
    h1mat=h1matX-h1matY;
    h2mat=h2matX-h2matY;
    params.DeltaTmat=reshape(arrayfun(@(h1,h2) DeltaT(h1,h2,params),h1mat(:),h2mat(:)),length(h1index),length(h1index));
    params.DeltaTTmat=reshape(arrayfun(@(h1,h2) DeltaTT(h1,h2,params),h1mat(:),h2mat(:)),length(h1index),length(h1index));
    params.Deltatmat=reshape(arrayfun(@(h1,h2) Deltal(h1,h2,-1,params),h1mat(:),h2mat(:)),length(h1index),length(h1index));
    params.Deltabmat=reshape(arrayfun(@(h1,h2) Deltal(h1,h2,1,params),h1mat(:),h2mat(:)),length(h1index),length(h1index));


    rotate=@(x) [cos(x) -sin(x);sin(x) cos(x)]; %rotate anticlockwise

    if params.nu==[2,1]
        params.q=[[0,0]];
        params.b=neighbor_index*[params.bM1;params.bM2];
        % line kappa_t-m-kappa_b-gamma
        m=(params.kb+params.kt)/2;
        kt_m_x=linspace(params.kt(1),m(1),40);
        kt_m_y=linspace(params.kt(2),m(2),40);
        m_kb_x=linspace(m(1),params.kb(1),40);
        m_kb_y=linspace(m(2),params.kb(2),40);
        kb_gamma_x=linspace(params.kb(1),0,40);
        kb_gamma_y=linspace(params.kb(2),0,40);


        kxlist=[kt_m_x,m_kb_x,kb_gamma_x];
        kylist=[kt_m_y,m_kb_y,kb_gamma_y];

        params.k=[kxlist',kylist'];

        segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
        params.klist=[0,cumsum(segment)];

        % params.k=[[0,0]];
    end