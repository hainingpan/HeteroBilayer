function [wbgrid,wtgrid]=w_rec(an,rx,ry,parameters)
%R: center of wannier state
%r: scalar|array real space position
%use rectangular (skew) grid
n=12;
state=1;
xrange=-n:n;
yrange=-n:n;
Nkx=length(xrange);
Nky=length(yrange);
% Nmax=parameters.Nmax;
bM1=parameters.bM1;
bM2=parameters.bM2;
%shift to diamond
% a1=-bM1/(2*n);
% a2=(bM1+bM2)/(2*n);
%shift to rectangular
a2=bM2/(2*n);
a1=(-2*bM1-bM2)/2/(2*n);

% [rx,ry]=meshgrid(linspace(-sqrt(3).aM,sqrt(3).aM,Nrx));
[Nrx,Nry]=size(rx);
% rectangular(skew grid)
kxmap=zeros(Nkx,Nky);
kymap=zeros(Nkx,Nky);
gauge=zeros(Nkx,Nky);
% expR=zeros(Nkx,Nky);
psibline=zeros(Nkx,Nky,Nrx*Nry);
psitline=zeros(Nkx,Nky,Nrx*Nry);

%Unit cell in momentum space
omega=abs(cross([bM1,0],[bM2,0]));
omega=omega(3);

%Unit cell in real space
V=abs(cross([parameters.aM1,0],[parameters.aM2,0]));
V=V(3);

% rectangular(skew) grid, for diamond mesh and rectangular mesh
for xindex=1:Nkx
    kx=xrange(xindex);
    for yindex=1:Nky
        ky=yrange(yindex);
        k=kx*a1+ky*a2;
        [~,vec]=energyTMD(k(1),k(2),parameters);
        psib0=sum(vec(1:end/2,state)); %bloch wf at r=(0,0) on the bottom layer
        gauge(xindex,yindex)=conj(abs(psib0)/psib0);
        [ubgrid,utgrid]=u(vec(:,state),rx,ry,parameters);
%         expR(xindex,yindex)=exp(-1i*dot(k,R));
        psibgrid=ubgrid.*exp(1i*(k(1)*rx+k(2)*ry))/sqrt(V);
        psitgrid=utgrid.*exp(1i*(k(1)*rx+k(2)*ry))/sqrt(V);
        psibline(xindex,yindex,:)=psibgrid(:);
        psitline(xindex,yindex,:)=psitgrid(:);        
        kxmap(xindex,yindex)=k(1);
        kymap(xindex,yindex)=k(2);
    end
end

for i=1:length(an)
    R=an{i}(1)*parameters.aM1+an{i}(2)*parameters.aM2;
    expR=exp(-1i*(kxmap*R(1)+kymap*R(2)));
    expRgauge=(gauge.*expR);
     
    %use integral, rectangular(skew) grid
    Fb=expRgauge.*psibline;
    Ft=expRgauge.*psitline;
    wbline=trapz(kxmap(:,1),trapz(kymap(1,:),Fb,2))/omega;
    wtline=trapz(kxmap(:,1),trapz(kymap(1,:),Ft,2))/omega;

    %use summation, rectangular(skew) grid,
    % psib=reshape(psibline,[Nkx*Nky,Nrx*Nry]);
    % psit=reshape(psitline,[Nkx*Nky,Nrx*Nry]);
    % wbline=expRgauge(:).'*psib/(Nkx*Nky);
    % wtline=expRgauge(:).'*psit/(Nkx*Nky);

    wbgrid(:,:,i)=reshape(wbline,[Nrx,Nry]);
    wtgrid(:,:,i)=reshape(wtline,[Nrx,Nry]);
end

end


