U = dlmread('U.dat');
th = linspace(0,2*pi,16);
r = linspace(0,1,11);
[TH,R] = meshgrid(th,r);
[X,Y] = pol2cart(TH,R);

figure
surf(X,Y,U)
xlabel('x')
ylabel('y')

figure
h = polar([0 2*pi], [0 1]);
delete(h)
hold on
contour(X,Y,U,30)