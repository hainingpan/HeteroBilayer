function neighborlist=generate_shell(nshell)
% For the usage of bM1, and bM2, which are 120 separated

counter=1;

for yindex=-nshell:nshell
    for xindex=max(-nshell,-nshell+yindex):min(nshell+yindex,nshell)
        x(counter)=xindex;
        y(counter)=yindex;
%         r(counter)=norm(xindex*a1+yindex*a2)^2;
        counter=counter+1;
    end
end

neighborlist=[x',y'];
% neighborlist=[[0,0];[0,1];[-1,0];[-1,-1];[0,-1];[1,0];[1,1]];
end
