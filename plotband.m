[kxq,kyq]=meshgrid(linspace(min(params_mK.k(:,1)),max(params_mK.k(:,1)),100),linspace(min(params_mK.k(:,2)),max(params_mK.k(:,2)),100));
energyall_sort=sort(energyall(:));
Nk=size(params_mK.k,1);
mu=energyall_sort(Nk*params_mK.nu(1)/(2*params_mK.nu(2)));
figure;
hold on;
surf(2*kxq/norm(params_mK.bm1),2*kyq/norm(params_mK.bm2),1000*mu*ones(100),'edgecolor','none','FaceAlpha',0.2);
% for i=1:size(energyall,2)
for i=1:3
    vq=griddata(params_mK.k(:,1),params_mK.k(:,2),1000*energyall(:,i),kxq,kyq);
    mesh(kxq/norm(params_mK.bm1),kyq/norm(params_mK.bm1),vq);
end
xlabel('k_x/|b_M|');
ylabel('k_y/|b_M|');
zlabel('E (meV)');
% title(sprintf("%d bands\ngap: %0.8f meV",size(energyall,2),1000*gap));
xlim([min(kxq(:)),max(kxq(:))]/norm(params_mK.bm1));
ylim([min(kyq(:)),max(kyq(:))]/norm(params_mK.bm1));

view([-1,-1,.1]);
