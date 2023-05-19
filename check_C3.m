function dH=check_C3(kx,ky,parameters)
[~,~,H]=energyTMD(kx,ky,parameters);
H11=H(1:end/2,1:end/2);
H12=H(1:end/2,end/2+1:end);
H21=H(end/2+1:end,1:end/2);
H22=H(end/2+1:end,end/2+1:end);

H11_0=H11;
H22_0=H22;
H12_0=H12;
H21_0=H21;

R=[0,-1;1,-1];
basis0=[parameters.h1index';parameters.h2index'];
basis1=R*basis0;
% basis0_q=basis0;
% basis1_q=basis1;
for index1=1:size(basis1,2)
    index0=1;
    while ~((basis1(1,index1)==basis0(1,index0)) && (basis1(2,index1)==basis0(2,index0)))
        index0=index0+1;
    end
    disp(index0)

    H11([index1,index0],:)=H11([index0,index1],:)
    H11(:,[index1,index0])=H11(:,[index0,index1],:)
    
    H22([index1,index0],:)=H22([index0,index1],:)
    H22(:,[index1,index0])=H22(:,[index0,index1],:)

    H12([index1,index0],:)=H12([index0,index1],:)
    H12(:,[index1,index0])=H12(:,[index0,index1],:)

    H21([index1,index0],:)=H21([index0,index1],:)
    H21(:,[index1,index0])=H21(:,[index0,index1],:)

H1=[H11,H12;H21,H22];
dH=H1-H;
end
    

        
