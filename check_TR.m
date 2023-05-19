function check_TR(H,params)
Nb=size(params.b,1);
% for k_index=1:length(params.k)
%     assert(norm(params.k_index(k_index,:)+params.k_index(end-k_index+1,:))<1e-10,sprintf('large error k at :%d',k_index))
% end

% for k_index=1:length(params.b)
%     assert(norm(params.b_index(k_index,:)+params.b_index(end-k_index+1,:))<1e-10,sprintf('large error b at :%d',k_index))
% end

for k_index=1:length(params.k)
    H_p=H(1:end/2,1:end/2,k_index);
    H_m=H(end/2+1:end,end/2+1:end,end-k_index+1);
    perm=[Nb:-1:1,2*Nb:-1:Nb+1];
    H_m_perm=H_m(perm,perm);
    dH=H_p-conj(H_m_perm);
    assert(sum(abs(dH(:)))<1e-10,sprintf('large error H at %d',k_index));
end

