function re=hubbardU_fft_2(wbgrid1,wtgrid1,neighborlist,rx,ry,parameters)
[Nx,Ny]=size(rx);
alpha=0.00729735; %e^2/(4*pi*epsilon0); c.f. Quicks Notes on oneNotes

Lx=rx(1,end)-rx(1,1);
Ly=ry(end,1)-ry(1,1);

kxlist=2*pi/Lx*(-floor(Nx/2):floor((Nx-1)/2));
kylist=2*pi/Ly*(-floor(Ny/2):floor((Ny-1)/2));

% wb1k=fftshift(fft2(abs(wbgrid1).^2))*Lx*Ly/((Nx-1)*(Ny-1));
% wt1k=fftshift(fft2(abs(wtgrid1).^2))*Lx*Ly/((Nx-1)*(Ny-1));
% 
% wb2k=reshape(wb1k(end:-1:1),Nx,Ny);
% wt2k=reshape(wt1k(end:-1:1),Nx,Ny);

Mq=fftshift(fft2(abs(wbgrid1).^2)+fft2(abs(wtgrid1).^2))*Lx*Ly/((Nx-1)*(Nx-1));


[kxmap,kymap]=meshgrid(kxlist,kylist);
% w2=reshape(wb1k.*wb2k+wt1k.*wt2k,Nx,Ny);
w2=abs(Mq).^2;

re=zeros(1,length(neighborlist));
F=griddedInterpolant(kxmap',kymap',w2);
parfor i=1:length(neighborlist)
    n=neighborlist{i}(1)*parameters.aM1+neighborlist{i}(2)*parameters.aM2;
    func=@(x,y) alpha/(2*pi)*1./(sqrt(x.^2+y.^2)).*(F(x',y')').*tanh(parameters.d*sqrt(x.^2+y.^2)).*exp(-1i*(x*n(1)+y*n(2)));
    re(i)=quad2d(func,kxlist(1),kxlist(end),kylist(1),kylist(end),'FailurePlot',false,'MaxFunEvals',5000);
end

end