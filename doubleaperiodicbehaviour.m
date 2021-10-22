clear
[t,y] = ode113(@vdp1,[0 50],[pi/2; pi/2; 0; 0]);
l = 1;
x1 = (l/2)*sin(y(:,1));
y1 = (-l/2)*cos(y(:,1));
x2 = (l/2)*(sin(y(:,1)) + sin(y(:,2))); y2 = (l/2)*(-cos(y(:,1)) - cos(y(:,2)));
plot (t, y(:,2)) 
xlabel('time') 
ylabel('theta 2')

function dydt = vdp1(t,y) 
g = 9.8;
l = 1;
m1 = 1;
m2 = 1;
dydt = [y(3); y(4); (((-2*g*(m1*sin(y(1)) + m2*sin(y(1))) + 2*m2*cos(y(1)-y(2))*g*sin(y(2)))/l) - m2*cos(y(1)-y(2))*y(3)*y(3)*sin(y(1)-y(2)) - m2*y(4)*y(4)*sin(y(1)-y(2)))/((m1+m2)-(m2*cos(y(1)-y(2))*cos(y(1)-y(2)))); ((- 2*g*sin(y(2)))/l + y(3)*y(3)*sin(y(1)-y(2)) + (2*g*cos(y(1)-y(2))*(m1*sin(y(1)) + m2*sin(y(1))) + m2*y(4)*y(4)*cos(y(1)-y(2))*sin(y(1)-y(2)))/(m1+m2))/(1-(m2*cos(y(1)-y(2))*cos(y(1)-y(2)))/(m1+m2))]
end
