[X,Y] = meshgrid(0:6,0:6);
U = 0.25*X;
V = 0.5*Y;
figure;quiver(X(:),Y(:),U(:),V(:),0)