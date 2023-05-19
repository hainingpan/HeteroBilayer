function U=U_calc_func(k,params)
% neighborlist=generate_shell(201);

xrange=-3*params.aM:params.aM/20:3*params.aM;
yrange=-3*params.aM:params.aM/20:3*params.aM;
[rx,ry]=meshgrid(xrange,yrange);

% neighborlist2=cellfun(@(x)x{1},neighborlist,'UniformOutput',false);
neighborlist2={[0,0],[0,1],[0,2]};
[wbgrid,wtgrid]=w_rec(neighborlist2(1),rx,ry,params);
% neighbordist=cellfun(@(x)norm(x*[parameters.aM1;parameters.aM2]),neighborlist2);



Uint=hubbardU_fft_2(wbgrid,wtgrid,[neighborlist2],rx,ry,params);
% for i=1:k+1
%     for j=1:length(neighborlist{i})
%         U{i}(j)=Utot(i);
%     end
% end
U=Uint;
end
