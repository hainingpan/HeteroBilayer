function [gap_list,tot_list,chern_p_list]=extract_title(nu,Nmax,epsilon,vz_t_list)
    gap_list=0*vz_t_list;
    tot_list=0*vz_t_list;
    for vz_t_index=1:length(vz_t_list)
       vz_t=vz_t_list(vz_t_index);
       fn=sprintf('D:\\CMTC\\Rp_HeteroBilayer\\3\\nu_%d,%d_Nmax%d_w%.1f_Nk_%d_Vzt_%.1f_ep%.1f_band.fig',nu,Nmax,20,15,vz_t,epsilon);
       fig=openfig(fn);
       h=gca;
       str1=sscanf(h.Title.String{1},'Gap: %e (meV) E: %e (meV)');
       gap=str1(1);
       tot=str1(2);
       gap_list(vz_t_index)=gap;
       tot_list(vz_t_index)=tot;
       str2=sscanf(h.Title.String{2},'(%d)+K:{%f} -K:{%f}');
       chern_p=str2(2);
       chern_p_list(vz_t_index)=chern_p;
    end
%     close all