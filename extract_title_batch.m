vz_t_list=-30:0;

for nu=[1]
    for ep=10:0.5:20
        for w=[12]
            for Vb=[7]
                [gap_list,tot_list,chern_p_list]=extract_title([nu,1],2,ep,w,Vb,5,vz_t_list);
                fn=sprintf('phase_nu%d,1_ep%.1f_w%.1f_Vb%.1f_d5.0_Nk15_gapfromband.mat',nu,ep,w,Vb);
                save(fn,'chern_p_list','gap_list','tot_list','vz_t_list');
                close all
            end
        end
    end
end
