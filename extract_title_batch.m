vz_t_list=-50:0;

for nu=[1,2]
    for ep=[10,15,20]
        for w=[4,12,20]
            for Vb=[5,10,15]
                [gap_list,tot_list,chern_p_list]=extract_title([nu,1],2,Vb,vz_t_list);
                fn=sprintf('phase_nu%d,1_ep%d_w%d_Vb%d_Nk_15_clean.mat',nu,ep,w,Vb);
                save(fn,'chern_p_list','gap_list','tot_list','vz_t_list');
                close all
            end
        end
    end
end