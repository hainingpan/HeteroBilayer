function [gap_indirect,gap_direct]=find_band_gap(params)
    en1_list=zeros(length(params.K),1);
    en2_list=zeros(length(params.K),1);
    parfor k_index=1:size(params.K,1)
        [en,~]=energyTMD(params.K(k_index,1),params.K(k_index,2),params);
        en1_list(k_index)=en(1);
        en2_list(k_index)=en(2);
    end
    gap_indirect=min(en1_list-en2_list);
    gap_direct=min(en1_list)-max(en2_list);
end