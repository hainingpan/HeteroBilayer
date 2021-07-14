n=10;
[xindexgrid,yindexgrid]=meshgrid(0:n-1,0:n-1);
xindexgrid/n.*params_mK.aM1+yindexgrid/n.*params_mK.aM2;