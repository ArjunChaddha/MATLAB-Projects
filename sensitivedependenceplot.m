clear
[t,y] = ode113(@vdp1,[0 50],[pi/2; 0]);
l = 1;
x1 = (l/2)*sin(y(:,1)); y1 = (-l/2)*cos(y(:,1));
plot (t, y(:,1)) xlabel('time') ylabel('theta 1')
function dydt = vdp1(t,y) g = 9.8;
l = 1;
m1 = 1;
m2 = 1;
dydt = [y(2); -(g/l)*(sin(y(1)))];
end
