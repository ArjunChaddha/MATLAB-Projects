clear
opts = odeset('RelTol',2.2205e-14,'AbsTol',2.2205e-14);
[t,y] = ode113(@vdp1,[0 40],[pi/2; pi/2; 0; 0], opts);

l = 1;
x1 = (l/2)*sin(y(:,1));
y1 = (-l/2)*cos(y(:,1));
x2 = (l/2)*(sin(y(:,1)) + sin(y(:,2)));
y2 = (l/2)*(-cos(y(:,1)) - cos(y(:,2)));

for k = 1:length(x1)
    clf
    h = animatedline('Marker', 'o');
    axis([-1.5,1.5,-1.5,1.5])
    axis square
    addpoints(h, 0, 0);
    addpoints(h,x1(k),y1(k));
    addpoints(h,x2(k),y2(k));
    drawnow
end

function dydt = vdp1(t,y)
    g = 9.8;
    l = 1;
    m1 = 1;
    m2 = 1;
    dydt = [y(3); y(4); (((-2*g*(m1*sin(y(1)) + m2*sin(y(1))) + 2*m2*cos(y(1)-y(2))*g*sin(y(2)))/l) - m2*cos(y(1)-y(2))*y(3)*y(3)*sin(y(1)-y(2)) - m2*y(4)*y(4)*sin(y(1)-y(2)))/((m1+m2)-(m2*cos(y(1)-y(2))*cos(y(1)-y(2)))); ((-2*g*sin(y(2)))/l + y(3)*y(3)*sin(y(1)-y(2)) + (2*g*cos(y(1)-y(2))*(m1*sin(y(1)) + m2*sin(y(1))) + m2*y(4)*y(4)*cos(y(1)-y(2))*sin(y(1)-y(2)))/(m1+m2))/(1-(m2*cos(y(1)-y(2))*cos(y(1)-y(2)))/(m1+m2))]
end