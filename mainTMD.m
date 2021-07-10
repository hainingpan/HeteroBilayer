function parameters=mainTMD(varargin)
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
    addParameter(p,'valley',-1);% -1: -K; +1; +k
    addParameter(p,'nu',[1,1]);%filling factor number(#1) per site(#2); no spin degeneracy considered
    addParameter(p,'hole',1); %1: hole-like energy band ; -1: particle-like energy band
    
    parse(p,varargin{:});
    parameters=struct('a_b',p.Results.a_b,'a_t',p.Results.a_t,'theta',p.Results.theta,'m_b',p.Results.m_b*0.511e6,'m_t',p.Results.m_t*0.511e6,'V_b',p.Results.V_b*1e-3,'V_t',p.Results.V_t*1e-3,'psi_b',p.Results.psi_b/360*2*pi,'psi_t',p.Results.psi_t/360*2*pi,'w',p.Results.w*1e-3,'Vz_b',p.Results.Vz_b*1e-3,'Vz_t',p.Results.Vz_t*1e-3,'Nmax',p.Results.Nmax,'valley',p.Results.valley,'nu',p.Results.nu,'hole',p.Results.hole);

    delta=(parameters.a_b-parameters.a_t)/parameters.a_t;
    parameters.aM=parameters.a_b/sqrt(delta^2+parameters.theta^2);
    parameters.aM1=parameters.aM*[1,0];
    parameters.aM2=parameters.aM*[1/2,sqrt(3)/2];
    parameters.g1=[0,4*pi/(sqrt(3)*parameters.aM)];
    parameters.g2=[-2*pi/parameters.aM,(2*pi)/(sqrt(3)*parameters.aM)];
    parameters.g3=[-2*pi/parameters.aM,-(2*pi)/(sqrt(3)*parameters.aM)];
    parameters.g4=[0,-4*pi/(sqrt(3)*parameters.aM)];
    parameters.g5=[2*pi/parameters.aM,-(2*pi)/(sqrt(3)*parameters.aM)];
    parameters.g6=[2*pi/parameters.aM,(2*pi)/(sqrt(3)*parameters.aM)];
    parameters.bM1=parameters.g5;
    parameters.bM2=parameters.g1;
%     parameters.kb=-4*pi/(3*parameters.aM)*[-1/2,sqrt(3)/2];
%     parameters.kt=-4*pi/(3*parameters.aM)*[1/2,sqrt(3)/2];
    parameters.kb=4*pi/(3*parameters.aM)*[-1/2,sqrt(3)/2];
    parameters.kt=4*pi/(3*parameters.aM)*[1/2,sqrt(3)/2];

    Nrange=-parameters.Nmax:parameters.Nmax;
    [h1index,h2index]=meshgrid(Nrange,Nrange);
    parameters.h1index=h1index;
    parameters.h2index=h2index;
    [h1matX,h1matY]=meshgrid(h1index(:));
    [h2matX,h2matY]=meshgrid(h2index(:));
    h1mat=h1matX-h1matY;
    h2mat=h2matX-h2matY;
    parameters.DeltaTmat=reshape(arrayfun(@(h1,h2) DeltaT(h1,h2,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);
    parameters.DeltaTTmat=reshape(arrayfun(@(h1,h2) DeltaTT(h1,h2,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);
    parameters.Deltatmat=reshape(arrayfun(@(h1,h2) Deltal(h1,h2,-1,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);
    parameters.Deltabmat=reshape(arrayfun(@(h1,h2) Deltal(h1,h2,1,parameters),h1mat(:),h2mat(:)),(2*parameters.Nmax+1)^2,(2*parameters.Nmax+1)^2);


