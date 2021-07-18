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
    addParameter(p,'d',5e-9*5.076e6); % gate to sample distance
    addParameter(p,'epsilon',1); % gate to sample distance
    addParameter(p,'Ez',0); % Additional potential to -K
    
    parse(p,varargin{:});
    params=struct('a_b',p.Results.a_b,'a_t',p.Results.a_t,'theta',p.Results.theta*pi/180,'m_b',p.Results.m_b*0.511e6,'m_t',p.Results.m_t*0.511e6,'V_b',p.Results.V_b*1e-3,'V_t',p.Results.V_t*1e-3,'psi_b',p.Results.psi_b/360*2*pi,'psi_t',p.Results.psi_t/360*2*pi,'w',p.Results.w*1e-3,'Vz_b',p.Results.Vz_b*1e-3,'Vz_t',p.Results.Vz_t*1e-3,'Nmax',p.Results.Nmax,'omega',p.Results.omega,'valley',p.Results.valley,'nu',p.Results.nu,'hole',p.Results.hole,'n',p.Results.n,'d',p.Results.d,'epsilon',p.Results.epsilon);

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
    %     params.shift=4*pi/(3*params.aM)*[1/2,sqrt(3)/2];
    % else
    %     params.shift=4*pi/(3*params.aM)*[-1/2,sqrt(3)/2];
    % end
    % params.kb=params.kb-params.shift;
    % params.kt=params.kt-params.shift;

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
    % params.DeltaTmat=reshape(arrayfun(@(h1,h2) DeltaT(h1,h2,params),h1mat(:),h2mat(:)),length(h1index),length(h1index));
    params.DeltaTmat=reshape(DeltaT(h1mat(:),h2mat(:),params),length(h1index),length(h1index));
    % params.DeltaTTmat=reshape(arrayfun(@(h1,h2) DeltaTT(h1,h2,params),h1mat(:),h2mat(:)),length(h1index),length(h1index));
    params.DeltaTTmat=reshape(DeltaTT(h1mat(:),h2mat(:),params),length(h1index),length(h1index));
    % params.Deltatmat=reshape(arrayfun(@(h1,h2) Deltal(h1,h2,-1,params),h1mat(:),h2mat(:)),length(h1index),length(h1index));
    params.Deltatmat=reshape(Deltal(h1mat(:),h2mat(:),-1,params),length(h1index),length(h1index));
    % params.Deltabmat=reshape(arrayfun(@(h1,h2) Deltal(h1,h2,1,params),h1mat(:),h2mat(:)),length(h1index),length(h1index));
    params.Deltabmat=reshape(Deltal(h1mat(:),h2mat(:),1,params),length(h1index),length(h1index));
    params.area=sqrt(3)/2*params.aM^2;

    if mod(params.n,3)~=0
        warning('n={%d} is not multiple of 3, which does not go through K point',params.n);
    end

    rotate=@(x) [cos(x) -sin(x);sin(x) cos(x)]; %rotate anticlockwise

    if params.nu==[1,1]
        params.q_index=[[0,0]];
        params.q=params.q_index*[params.bM1;params.bM2];

        params.b_index=neighbor_index;
        params.b=params.b_index*[params.bM1;params.bM2];

        params.bm1=params.bM1;
        params.bm2=params.bM2;
        
        [ux,uy]=ndgrid(1:params.n,1:params.n);
        params.k_index=[(2*ux(:)-2)/(2*params.n),(2*uy(:)-2)/(2*params.n)];
        params.k=params.k_index*[params.bM1;params.bM2];
    end

    if params.nu==[2,1]
        params.q_index=[[0,0]];
        params.q=params.q_index*[params.bM1;params.bM2];

        params.b_index=neighbor_index;
        params.b=params.b_index*[params.bM1;params.bM2];

        params.bm1=params.bM1;
        params.bm2=params.bM2;
        
        [ux,uy]=ndgrid(1:params.n,1:params.n);
        % params.k_index=[(2*ux(:)-params.n-1)/(2*params.n),(2*uy(:)-params.n-1)/(2*params.n)];
        params.k_index=[(2*ux(:)-params.n)/(2*params.n),(2*uy(:)-params.n)/(2*params.n)];
        % params.k_index=[(2*ux(:)-2)/(2*params.n),(2*uy(:)-2)/(2*params.n)];
        % params.k_index=[(2*ux(:)-7)/(2*params.n),(2*uy(:)-7)/(2*params.n)];
        params.k=params.k_index*[params.bM1;params.bM2];
        


        % line kappa_t-m-kappa_b-gamma
        % m=(params.kb+params.kt)/2;
        % kt_m_x=linspace(params.kt(1),m(1),40);
        % kt_m_y=linspace(params.kt(2),m(2),40);
        % m_kb_x=linspace(m(1),params.kb(1),40);
        % m_kb_y=linspace(m(2),params.kb(2),40);
        % kb_gamma_x=linspace(params.kb(1),0,40);
        % kb_gamma_y=linspace(params.kb(2),0,40);


        % kxlist=[kt_m_x,m_kb_x,kb_gamma_x];
        % kylist=[kt_m_y,m_kb_y,kb_gamma_y];

        % params.k=[kxlist',kylist'];

        % segment=sqrt(diff(kxlist).^2+diff(kylist).^2);
        % params.klist=[0,cumsum(segment)];
    end



    % am1=am1index*[params.aM1;params.aM2];
    % am2=am2index*[params.aM1;params.aM2];
    % params.bm1=2*pi/abs(length(ailist)*params.aM^2*sqrt(3)/2)*am1*rotate(-pi/2);
    % params.bm2=2*pi/abs(length(ailist)*params.aM^2*sqrt(3)/2)*am2*rotate(pi/2);
    % qindex=[params.bm1;params.bm2]/[params.bM1;params.bM2];
    % manhattandist=1;
    % Qlist={[0,0]};
    % done=0;
    % while 1  %Very ugly code
    %     for i=-manhattandist:manhattandist
    %         if abs(i)==manhattandist
    %             jlist=[0];
    %         else
    %             jlist=[-(manhattandist-abs(i)),(manhattandist-abs(i))];
    %         end
    %         for j=jlist
    %         qtmp=[i,j]*qindex;
    %         if all(qtmp>=-1e-10) & all(qtmp<1-1e-10)
    %             Qlist=[Qlist,qtmp];
    %         end
    %         if length(Qlist)==length(ailist)
    %             done=1;
    %             break;
    %         end
    %         end
    %     end
    %     if done==1
    %         break;
    %     end
    %     manhattandist=manhattandist+1;
    % end
    % params.Qindex=Qlist;
    % params.Q=cellfun(@(x) x(1)*params.bM1+x(2)*params.bM2,Qlist,'UniformOutput',0);

    if params.nu~=0
        [b_a_x,b_b_x]=meshgrid(params.b_index(:,1),params.b_index(:,1));
        [b_a_y,b_b_y]=meshgrid(params.b_index(:,2),params.b_index(:,2));
        params.valley=1;
        params.Delta_T_p=reshape(DeltaT(b_a_x-b_b_x,b_a_y-b_b_y,params),size(params.b_index,1),size(params.b_index,1));
        params.Delta_TT_p=reshape(DeltaTT(b_a_x-b_b_x,b_a_y-b_b_y,params),size(params.b_index,1),size(params.b_index,1));
        params.Delta_b_p=reshape(Deltal(b_a_x-b_b_x,b_a_y-b_b_y,1,params),size(params.b_index,1),size(params.b_index,1));
        params.Delta_t_p=reshape(Deltal(b_a_x-b_b_x,b_a_y-b_b_y,-1,params),size(params.b_index,1),size(params.b_index,1));
        params.valley=-1;
        params.Delta_T_m=reshape(DeltaT(b_a_x-b_b_x,b_a_y-b_b_y,params),size(params.b_index,1),size(params.b_index,1));
        params.Delta_TT_m=reshape(DeltaTT(b_a_x-b_b_x,b_a_y-b_b_y,params),size(params.b_index,1),size(params.b_index,1));
        params.Delta_b_m=reshape(Deltal(b_a_x-b_b_x,b_a_y-b_b_y,1,params),size(params.b_index,1),size(params.b_index,1));
        params.Delta_t_m=reshape(Deltal(b_a_x-b_b_x,b_a_y-b_b_y,-1,params),size(params.b_index,1),size(params.b_index,1));


        % For the detailed def, check oneNotes
        [q_a_x,q_b_x,q_g_x,q_d_x,b_a_x,b_b_x,b_g_x,b_d_x]=ndgrid(params.q_index(:,1),params.q_index(:,1),params.q_index(:,1),params.q_index(:,1),params.b_index(:,1),params.b_index(:,1),params.b_index(:,1),params.b_index(:,1));
        [q_a_y,q_b_y,q_g_y,q_d_y,b_a_y,b_b_y,b_g_y,b_d_y]=ndgrid(params.q_index(:,2),params.q_index(:,2),params.q_index(:,2),params.q_index(:,2),params.b_index(:,2),params.b_index(:,2),params.b_index(:,2),params.b_index(:,2));
        delta_x=-q_a_x-b_a_x+q_b_x+b_b_x+q_g_x+b_g_x-q_d_x-b_d_x;
        delta_y=-q_a_y-b_a_y+q_b_y+b_b_y+q_g_y+b_g_y-q_d_y-b_d_y;
        delta_tensor=abs(delta_x)<1e-10 & abs(delta_y)<1e-10;
        params.delta_tensor1=int8(delta_tensor); %-q_a-b_a+q_b+b_b+q_g+b_g-q_d-b_d, index: q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d

        delta_x=-q_a_x-b_a_x+q_b_x+b_b_x-q_g_x-b_g_x+q_d_x+b_d_x;
        delta_y=-q_a_y-b_a_y+q_b_y+b_b_y-q_g_y-b_g_y+q_d_y+b_d_y;
        delta_tensor=abs(delta_x)<1e-10 & abs(delta_y)<1e-10;
        params.delta_tensor2=int8(delta_tensor); %-q_a-b_a+q_b+b_b-q_g-b_g+q_d+b_d, index: q_a,q_b,q_g,q_d,b_a,b_b,b_g,b_d

        alpha=0.00729735; % eV*nm

        [q_g_x,q_d_x,b_g_x,b_d_x]=ndgrid(params.q(:,1),params.q(:,1),params.b(:,1),params.b(:,1));
        [q_g_y,q_d_y,b_g_y,b_d_y]=ndgrid(params.q(:,2),params.q(:,2),params.b(:,2),params.b(:,2));
        q_abs=sqrt((b_g_x+q_g_x-b_d_x-q_d_x).^2+(b_g_y+q_g_y-b_d_y-q_d_y).^2);
        qd=q_abs*params.d;
        % params.V1=alpha/2*sinc(1i*qd/pi)./cos(1i*qd)*params.d; % a work-around to avoid 0/0, index: q_g,q_d,b_g,b_d
        %         params.V1(isnan(params.V1))=1/2*params.d;
        params.V1=alpha/2*tanh(qd+1e-18)./(qd+1e-18)*params.d;
        assert(sum(isnan(params.V1(:)))==0,'V1 has NaN');

        [k_a_x,k_b_x,q_a_x,q_d_x,b_a_x,b_d_x]=ndgrid(params.k(:,1),params.k(:,1),params.q(:,1),params.q(:,1),params.b(:,1),params.b(:,1));
        [k_a_y,k_b_y,q_a_y,q_d_y,b_a_y,b_d_y]=ndgrid(params.k(:,2),params.k(:,2),params.q(:,2),params.q(:,2),params.b(:,2),params.b(:,2));
        q_abs=sqrt((k_a_x+b_d_x+q_d_x-k_b_x-b_a_x-q_a_x).^2+(k_a_y+b_d_y+q_d_y-k_b_y-b_a_y-q_a_y).^2);
        qd=q_abs*params.d;
        % params.V2=alpha/2*sinc(1i*qd/pi)./cos(1i*qd)*params.d; % a work-around to avoid 0/0, index: k_a,k_b,q_a,q_d,b_a,b_d
        params.V2=alpha/2*tanh(qd+1e-18)./(qd+1e-18)*params.d;
        %params.V2(isnan(params.V2))=1/2*params.d;
        assert(sum(isnan(params.V2(:)))==0,'V2 has NaN');

        qindex=[params.bm1;params.bm2]/[params.bM1;params.bM2];


        Qlistmat=params.q;
        Qshift1=mod(Qlistmat+qindex(1,:),1);
        Qshift2=mod(Qlistmat+qindex(2,:),1);

        Qshift1=mod(Qshift1,1);
        Qshift2=mod(Qshift2,1);

        perm1=zeros(size(Qshift1,1));
        for i=1:size(Qshift1,1)
            j=find(abs(sum(abs(Qlistmat(i,:)-Qshift1).^2,2))<1e-5);
            perm1(i,j)=1;
        end

        perm2=zeros(size(Qshift2,1));
        for i=1:size(Qshift2,1)
            j=find(abs(sum(abs(Qlistmat(i,:)-Qshift2).^2,2))<1e-5);
            perm2(i,j)=1;
        end 

        params.perm1=perm1;
        params.perm2=perm2;



    end