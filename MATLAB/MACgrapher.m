% navigate to save directory (containing u,v,p files) and run this code to
% plot the Navier-Stokes simulation

x0 = 0;
xn = 1;
y0 = 0;
yn = 1;

n = 64;
h = (xn-x0)/n;
x = x0+h/2:h:xn-h/2;
y = (y0+h/2:h:yn-h/2);


hfg = figure('Position', [100 100 900 750]);
% axis tight manual
% filename = 'testAnimated.gif';
for i = 1:1:400
Ax = load(sprintf("u/%d.txt",i-1));
Ay = load(sprintf("v/%d.txt",i-1));
Ap = load(sprintf("p/%d.txt",i-1));

Ux = reshape(Ax',n,n)';
Uy = reshape(Ay',n,n)';
Up = reshape(Ap',n,n)';
        
[X,Y] = meshgrid(x,y);

clf

shading interp;
imagesc(x,(x0+h/2):h:(xn-h/2),Up);
set(gca,'YDir','normal');
hold on;
plot(0.7+(0.15)*cos((2*pi/128)*(0:128)),0.5+(0.15)*sin((2*pi/128)*(0:128)),"k","LineWidth",2);
hh = streamslice(X,Y,Ux,Uy,3);
set(hh,'Color',[0 0 0]);


colorbar;
caxis([0 0.75]);
xlim([x0,xn])
ylim([y0,yn])

title(i-1)      
drawnow
% frame = getframe(hfg);
% im = frame2im(frame);
% [imind,cm] = rgb2ind(im,256);
% if i == 1
%   imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
% else
%   imwrite(imind,cm,filename,'gif','WriteMode','append');
% end
% pause on
% pause(1.0)
end