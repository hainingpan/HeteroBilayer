n=30;
aa=[1/3,1/3;-2/3,1/3]*[params.aM1/n;params.aM2/n];
rr=generate_shell(n)*aa;
rr=rr+2*params.aM1;
[s0,sx,sy,sz]=S_r(ave2,rr(:,1),rr(:,2),1,params);