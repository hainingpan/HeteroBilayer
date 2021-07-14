function ave=average(energyall,wfall,params)
    energyall_sort=sort(energyall(:));
    mu=energyall_sort(end*parameters.nu(1)/(2*parameters.nu(2)));
