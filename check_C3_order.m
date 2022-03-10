function check_C3_order(params)
    % The order parameter is defined as |<psi|C3|psi>|. 1 is C3 preserved.
    [energyall_p,energyall_m,wfall_p,wfall_m]=energyMF_bc(ave1,ave2,'C3',epoch,params)

