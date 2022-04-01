function [dos_list,en_list,nu_list,E_vanHove]=DOS(energy_list,params)
    en_list=linspace(min(energy_list),max(energy_list),200);
    eta=1e-4;
    dos_list=zeros(1,length(en_list));
    nu_list=zeros(1,length(en_list));
    for i=1:length(en_list)
        deltaf=1/pi*eta./((en_list(i)-energy_list).^2+eta^2);
        dos_list(i)=sum(deltaf(:));
        nu_list(i)=sum(energy_list>=en_list(i))/length(energy_list);
    end
    dos_list=dos_list/length(energy_list)/(sqrt(3)/2*(params.aM/5.076e-3)^2);
    [~,I]=max(dos_list);
    E_vanHove=en_list(I);
end