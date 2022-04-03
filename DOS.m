function [dos_list,en_list,nu_list,E_vanHove]=DOS(energy_list,params)
    % [dos_list,en_list,nu_list,E_vanHove]=
    % en_list=linspace(min(energy_list),max(energy_list),200);
    % eta=1e-4;
    % dos_list=zeros(1,length(en_list));
    % nu_list=zeros(1,length(en_list));
    % for i=1:length(en_list)
    %     deltaf=1/pi*eta./((en_list(i)-energy_list).^2+eta^2);
    %     dos_list(i)=sum(deltaf(:));
    %     nu_list(i)=sum(energy_list>=en_list(i))/length(energy_list);
    % end
    % dos_list=dos_list/length(energy_list)/(sqrt(3)/2*(params.aM/5.076e-3)^2);
    % [~,I]=max(dos_list);
    % E_vanHove=en_list(I);

    % Omega number of states
    energy_list=energy_list(~isnan(energy_list));
    energy_list=sort(energy_list);
    [dos_list,en_list]=ksdensity(energy_list,'bandwidth',0.1e-3);
    dos_list=dos_list/(sqrt(3)/2*(params.aM/5.076e-3)^2);
    [~,I]=max(dos_list);
    E_vanHove=en_list(I);
    
    Omega_list=zeros(length(en_list),1);
    for i=1:length(en_list)
        Omega_list(i)=sum(en_list(i)>=energy_list);
    end
    nu_list=Omega_list/length(energy_list);
 
end