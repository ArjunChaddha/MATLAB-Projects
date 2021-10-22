clear
theta1 = pi/2;
theta2 = pi/2;
T = 30;
numOfPts = 36;
 
TT = linspace (0, T, 10000);
 
opts = odeset('RelTol', 2.2205e-14, 'AbsTol', 2.2205e-14);
[t,y] = ode113(@vdp1,[TT],[theta1; theta2; 0; 0], opts);
 
angles = linspace(-pi, pi, numOfPts);
 
s = -13;
e = -4;
num = -1*(s-e-1);
xx = zeros(1, num) + 10;
radius = linspace(s,e, num);
radius = xx.^radius;
 
Ths = zeros(1, length(radius));
for r=1:length(radius)
    matrixTh1 = zeros(1,length(TT));
    matrixTh2 = zeros(1,length(TT));
 
    for i=1:length(angles)
        newtheta1 = theta1 + (radius(r))*cos(angles(i));
        newtheta2 = theta2 + (radius(r))*sin(angles(i));
        [t1, y1] = ode113(@vdp1, [TT], [newtheta1; newtheta2; 0; 0], opts);
 
 
        deltatheta1 = (y1(:,1)-y(:,1));
        deltatheta2 = (y1(:,2)-y(:,2));
        deltatheta1 = transpose(deltatheta1);
        deltatheta2 = transpose(deltatheta2);
        
        matrixTh1 = [matrixTh1; deltatheta1];
        matrixTh2 = [matrixTh2; deltatheta2];
    end
 
 
    maxS = zeros(1,length(TT));
    for k = 1:length(TT)
        dists = zeros(1, numOfPts);
        for i = 1:length(angles)
            distance = ((matrixTh1(i, k))^2 + (matrixTh2(i,k)^2))^0.5;
            dists(i) = distance;
        end
        maximum = max(dists);
        maxS(k) = maximum;
    end
 
    subtracted = (maxS - 10^-3);

    counter = 0;
    flag = 0;
    for t = 1:length(TT)
        if subtracted(t)<=0 & flag == 0
            counter = counter + 1;
        end
        if subtracted(t)>0
            flag = 1;
        end
    end
    
    TT(counter)
    radius(r)
    Ths(r) = TT(counter);     
        

end

plot(radius, Ths)
plot(log10(radius), Ths)
 
function dydt = vdp1(t,y)
    g = 9.8;
    l = 1;
    m1 = 1;
    m2 = 1;
    dydt = [y(3); y(4); (((-2*g*(m1*sin(y(1)) + m2*sin(y(1))) + 2*m2*cos(y(1)-y(2))*g*sin(y(2)))/l) - m2*cos(y(1)-y(2))*y(3)*y(3)*sin(y(1)-y(2)) - m2*y(4)*y(4)*sin(y(1)-y(2)))/((m1+m2)-(m2*cos(y(1)-y(2))*cos(y(1)-y(2)))); ((-2*g*sin(y(2)))/l + y(3)*y(3)*sin(y(1)-y(2)) + (2*g*cos(y(1)-y(2))*(m1*sin(y(1)) + m2*sin(y(1))) + m2*y(4)*y(4)*cos(y(1)-y(2))*sin(y(1)-y(2)))/(m1+m2))/(1-(m2*cos(y(1)-y(2))*cos(y(1)-y(2)))/(m1+m2))];
end