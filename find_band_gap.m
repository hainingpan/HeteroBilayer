function en_min=find_band_gap(params)
    en1_list=zeros(length(params.K),1);
    en2_list=zeros(length(params.K),1);
    parfor k_index=1:length(params.K)
        [en,~]=energyTMD(params.K(k_index,1),params.K(k_index,2),params);
        en1_list(k_index)=en(1);
        en2_list(k_index)=en(2);
    end
    % en_min=min(en1_list-en2_list);
    en_min=min(en1_list)-max(en2_list);
end