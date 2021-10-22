opts = odeset('RelTol',2.2205e-14,'AbsTol',2.2205e-14);
[t,y] = ode113(@vdp1,[0 40],[pi+0.5; 0], opts)

 
 
h = animatedline;
axis([-1.5,1.5,-1,0])
 
x = sin(y(:,1));
y = -cos(y(:,1));

for k = 1:length(x1)
    clf
    h = animatedline('Marker', 'o');
    axis([-1.5,1.5,-1.5,1.5])
    axis square
    addpoints(h, 0, 0);
    addpoints(h,x(k),y(k));
    drawnow
end
 
function dydt = vdp1(t,y)
    g = 9.8;
    l = 1;
    dydt = [y(2); (-9.8*sin(y(1)))/l];
end
 