function S_r_wf(rx,ry,wfall,params)
wf_b_p=wfall(:,1:2,1:end/4);
wf_t_p=wfall(:,1:2,end/4+1:end/2);
wf_b_m=wfall(:,1:2,end/2+1:end/4*3);
wf_t_m=wfall(:,1:2,end/4*3+1:end);

[q_x,b_x,r_x]=ndgrid(params.q(:,1),params.b(:,1),rx);
[q_y,b_y,r_y]=ndgrid(params.q(:,2),params.b(:,2),ry);
expm=exp(1i*((q_x+b_x).*r_x+(q_y+b_y).*r_y)); % q,b,r



