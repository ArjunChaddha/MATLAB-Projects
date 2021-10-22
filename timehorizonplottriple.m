clear
theta1 = pi/2;
theta2 = pi/2;
theta3 = pi/2;
T = 15;
TT = linspace (0, T, 5000);
numOfPts = 18;
opts = odeset('RelTol', 2.2205e-14, 'AbsTol', 2.2205e-14);
[t,y] = ode113(@vdp1,[TT],[theta1; 0; theta2; 0; theta3; 0], opts);
 
 
angles1 = linspace(-pi, pi, numOfPts);
angles2 = linspace(-pi, pi, numOfPts);
 
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
    matrixTh3 = zeros(1, length(TT));
    o=1
    for i=1:length(angles1)
        for j=1:length(angles2)
            newtheta1 = theta1 + (radius(r))*sin(angles1(i))*cos(angles2(j));
            newtheta2 = theta2 + (radius(r))*sin(angles1(i))*sin(angles2(j));
            newtheta3 = theta3 + (radius(r))*cos(angles1(i));
            [t1, y1] = ode113(@vdp1, [TT], [newtheta1; 0; newtheta2; 0; newtheta3; 0], opts);
   
            deltatheta1 = (y1(:,1)-y(:,1));
            deltatheta2 = (y1(:,3)-y(:,3));
            deltatheta3 = (y1(:,5)-y(:,5));
            deltatheta1 = transpose(deltatheta1);
            deltatheta2 = transpose(deltatheta2);
            deltatheta3 = transpose(deltatheta3);
 
            matrixTh1 = [matrixTh1; deltatheta1];
            matrixTh2 = [matrixTh2; deltatheta2];
            matrixTh3 = [matrixTh3; deltatheta3];
        end
    end
 

 
    maxS = zeros(1,length(TT));
    for k = 1:length(TT)
        dists = zeros(1, numOfPts);
        for i = 1:length(angles1)
            distance = ((matrixTh1(i, k))^2 + (matrixTh2(i,k)^2) + (matrixTh3(i,k)^2))^0.5;
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
        
 
    % plot3(matrixTh1(:,:), matrixTh2(:, :), TT)
    % xlabel('x(t)')
    % ylabel('y(t)')
    % zlabel('z(t)')
    % grid on
end
 
plot(log10(radius), Ths)
 
function dydt = vdp1(t,y)
    g = 9.8;
    l = 1;
    m1 = 1;
    m2 = 1;
    dydt = [y(2); -(588*cos(y(3))^2*cos(y(5))^2*sin(y(1)) + 441*cos(y(3))^2*sin(y(1))*sin(y(5))^2 - 441*cos(y(1))*cos(y(3))*cos(y(5))^2*sin(y(3)) - 147*cos(y(1))*cos(y(3))^2*cos(y(5))*sin(y(5)) - 588*cos(y(1))*cos(y(3))*sin(y(3))*sin(y(5))^2 - 15*cos(y(1))*cos(y(5))^2*sin(y(3))^3*y(4)^2 + 20*cos(y(3))^3*cos(y(5))^2*sin(y(1))*y(4)^2 - 5*cos(y(1))*cos(y(3))^2*sin(y(5))^3*y(6)^2 + 10*cos(y(3))^2*cos(y(5))^3*sin(y(1))*y(6)^2 - 20*cos(y(1))*sin(y(3))^3*sin(y(5))^2*y(4)^2 + 15*cos(y(3))^3*sin(y(1))*sin(y(5))^2*y(4)^2 - 10*cos(y(1))*sin(y(3))^2*sin(y(5))^3*y(6)^2 + 5*cos(y(5))^3*sin(y(1))*sin(y(3))^2*y(6)^2 + 147*cos(y(3))*cos(y(5))*sin(y(1))*sin(y(3))*sin(y(5)) - 5*cos(y(1))*cos(y(3))*cos(y(5))^3*sin(y(3))*y(6)^2 - 5*cos(y(1))*cos(y(3))^3*cos(y(5))*sin(y(5))*y(4)^2 + 5*cos(y(3))*sin(y(1))*sin(y(3))*sin(y(5))^3*y(6)^2 + 5*cos(y(5))*sin(y(1))*sin(y(3))^3*sin(y(5))*y(4)^2 + 15*cos(y(1))*cos(y(3))^2*cos(y(5))^2*sin(y(1))*y(2)^2 - 10*cos(y(1))^2*cos(y(3))*cos(y(5))^2*sin(y(3))*y(2)^2 - 15*cos(y(1))*cos(y(3))^2*cos(y(5))^2*sin(y(3))*y(4)^2 - 5*cos(y(1))^2*cos(y(3))^2*cos(y(5))*sin(y(5))*y(2)^2 - 5*cos(y(1))*cos(y(3))^2*cos(y(5))^2*sin(y(5))*y(6)^2 + 5*cos(y(1))*cos(y(3))^2*sin(y(1))*sin(y(5))^2*y(2)^2 - 5*cos(y(1))*cos(y(5))^2*sin(y(1))*sin(y(3))^2*y(2)^2 + 10*cos(y(3))*cos(y(5))^2*sin(y(1))^2*sin(y(3))*y(2)^2 - 10*cos(y(1))^2*cos(y(3))*sin(y(3))*sin(y(5))^2*y(2)^2 - 20*cos(y(1))*cos(y(3))^2*sin(y(3))*sin(y(5))^2*y(4)^2 + 20*cos(y(3))*cos(y(5))^2*sin(y(1))*sin(y(3))^2*y(4)^2 - 5*cos(y(1))^2*cos(y(5))*sin(y(3))^2*sin(y(5))*y(2)^2 + 5*cos(y(3))^2*cos(y(5))*sin(y(1))^2*sin(y(5))*y(2)^2 - 10*cos(y(1))*cos(y(5))^2*sin(y(3))^2*sin(y(5))*y(6)^2 + 10*cos(y(3))^2*cos(y(5))*sin(y(1))*sin(y(5))^2*y(6)^2 - 15*cos(y(1))*sin(y(1))*sin(y(3))^2*sin(y(5))^2*y(2)^2 + 10*cos(y(3))*sin(y(1))^2*sin(y(3))*sin(y(5))^2*y(2)^2 + 15*cos(y(3))*sin(y(1))*sin(y(3))^2*sin(y(5))^2*y(4)^2 + 5*cos(y(5))*sin(y(1))^2*sin(y(3))^2*sin(y(5))*y(2)^2 + 5*cos(y(5))*sin(y(1))*sin(y(3))^2*sin(y(5))^2*y(6)^2 - 5*cos(y(1))*cos(y(3))*cos(y(5))*sin(y(3))^2*sin(y(5))*y(4)^2 - 5*cos(y(1))*cos(y(3))*cos(y(5))*sin(y(3))*sin(y(5))^2*y(6)^2 + 5*cos(y(3))^2*cos(y(5))*sin(y(1))*sin(y(3))*sin(y(5))*y(4)^2 + 5*cos(y(3))*cos(y(5))^2*sin(y(1))*sin(y(3))*sin(y(5))*y(6)^2)/(5*(sin(y(1))^2*sin(y(3))^2*sin(y(5))^2 + cos(y(1))^2*cos(y(3))^2*cos(y(5))^2 + 2*cos(y(1))^2*cos(y(3))^2*sin(y(5))^2 + 3*cos(y(1))^2*cos(y(5))^2*sin(y(3))^2 + 4*cos(y(3))^2*cos(y(5))^2*sin(y(1))^2 + 4*cos(y(1))^2*sin(y(3))^2*sin(y(5))^2 + 3*cos(y(3))^2*sin(y(1))^2*sin(y(5))^2 + 2*cos(y(5))^2*sin(y(1))^2*sin(y(3))^2 - 4*cos(y(1))*cos(y(3))*cos(y(5))^2*sin(y(1))*sin(y(3)) - 2*cos(y(1))*cos(y(3))^2*cos(y(5))*sin(y(1))*sin(y(5)) - 4*cos(y(1))*cos(y(3))*sin(y(1))*sin(y(3))*sin(y(5))^2 - 2*cos(y(1))*cos(y(5))*sin(y(1))*sin(y(3))^2*sin(y(5)))); y(4); -(441*cos(y(1))^2*cos(y(5))^2*sin(y(3)) + 588*cos(y(1))^2*sin(y(3))*sin(y(5))^2 + 294*cos(y(5))^2*sin(y(1))^2*sin(y(3)) + 147*sin(y(1))^2*sin(y(3))*sin(y(5))^2 - 294*cos(y(1))*cos(y(3))*cos(y(5))^2*sin(y(1)) - 294*cos(y(1))*cos(y(3))*sin(y(1))*sin(y(5))^2 - 10*cos(y(3))*cos(y(5))^2*sin(y(1))^3*y(2)^2 + 10*cos(y(1))^3*cos(y(5))^2*sin(y(3))*y(2)^2 + 5*cos(y(1))^2*cos(y(5))^3*sin(y(3))*y(6)^2 - 10*cos(y(3))*sin(y(1))^3*sin(y(5))^2*y(2)^2 + 10*cos(y(1))^3*sin(y(3))*sin(y(5))^2*y(2)^2 - 5*cos(y(3))*sin(y(1))^2*sin(y(5))^3*y(6)^2 - 294*cos(y(1))*cos(y(5))*sin(y(1))*sin(y(3))*sin(y(5)) - 5*cos(y(1))*cos(y(3))*cos(y(5))^3*sin(y(1))*y(6)^2 + 5*cos(y(1))*sin(y(1))*sin(y(3))*sin(y(5))^3*y(6)^2 - 10*cos(y(1))^2*cos(y(3))*cos(y(5))^2*sin(y(1))*y(2)^2 - 10*cos(y(1))*cos(y(3))^2*cos(y(5))^2*sin(y(1))*y(4)^2 + 10*cos(y(1))^2*cos(y(3))*cos(y(5))^2*sin(y(3))*y(4)^2 + 10*cos(y(1))*cos(y(5))^2*sin(y(1))^2*sin(y(3))*y(2)^2 - 10*cos(y(1))^2*cos(y(3))*sin(y(1))*sin(y(5))^2*y(2)^2 - 10*cos(y(1))*cos(y(3))^2*sin(y(1))*sin(y(5))^2*y(4)^2 + 10*cos(y(1))*cos(y(5))^2*sin(y(1))*sin(y(3))^2*y(4)^2 - 10*cos(y(3))*cos(y(5))^2*sin(y(1))^2*sin(y(3))*y(4)^2 + 10*cos(y(1))^2*cos(y(3))*sin(y(3))*sin(y(5))^2*y(4)^2 - 5*cos(y(3))*cos(y(5))^2*sin(y(1))^2*sin(y(5))*y(6)^2 + 5*cos(y(1))^2*cos(y(5))*sin(y(3))*sin(y(5))^2*y(6)^2 + 10*cos(y(1))*sin(y(1))^2*sin(y(3))*sin(y(5))^2*y(2)^2 + 10*cos(y(1))*sin(y(1))*sin(y(3))^2*sin(y(5))^2*y(4)^2 - 10*cos(y(3))*sin(y(1))^2*sin(y(3))*sin(y(5))^2*y(4)^2 - 5*cos(y(1))*cos(y(3))*cos(y(5))*sin(y(1))*sin(y(5))^2*y(6)^2 + 5*cos(y(1))*cos(y(5))^2*sin(y(1))*sin(y(3))*sin(y(5))*y(6)^2)/(5*(sin(y(1))^2*sin(y(3))^2*sin(y(5))^2 + cos(y(1))^2*cos(y(3))^2*cos(y(5))^2 + 2*cos(y(1))^2*cos(y(3))^2*sin(y(5))^2 + 3*cos(y(1))^2*cos(y(5))^2*sin(y(3))^2 + 4*cos(y(3))^2*cos(y(5))^2*sin(y(1))^2 + 4*cos(y(1))^2*sin(y(3))^2*sin(y(5))^2 + 3*cos(y(3))^2*sin(y(1))^2*sin(y(5))^2 + 2*cos(y(5))^2*sin(y(1))^2*sin(y(3))^2 - 4*cos(y(1))*cos(y(3))*cos(y(5))^2*sin(y(1))*sin(y(3)) - 2*cos(y(1))*cos(y(3))^2*cos(y(5))*sin(y(1))*sin(y(5)) - 4*cos(y(1))*cos(y(3))*sin(y(1))*sin(y(3))*sin(y(5))^2 - 2*cos(y(1))*cos(y(5))*sin(y(1))*sin(y(3))^2*sin(y(5)))); y(6); -((cos(y(1))*sin(y(5)) - cos(y(5))*sin(y(1)))*(294*cos(y(1))*cos(y(3))^2 + 10*cos(y(1))*cos(y(3))^3*y(4)^2 + 10*sin(y(1))*sin(y(3))^3*y(4)^2 + 10*cos(y(1))^2*cos(y(3))^2*y(2)^2 + 10*cos(y(1))^2*sin(y(3))^2*y(2)^2 + 10*cos(y(3))^2*sin(y(1))^2*y(2)^2 + 10*sin(y(1))^2*sin(y(3))^2*y(2)^2 + 294*cos(y(3))*sin(y(1))*sin(y(3)) + 5*cos(y(1))*cos(y(3))^2*cos(y(5))*y(6)^2 + 10*cos(y(1))*cos(y(3))*sin(y(3))^2*y(4)^2 + 5*cos(y(1))*cos(y(5))*sin(y(3))^2*y(6)^2 + 10*cos(y(3))^2*sin(y(1))*sin(y(3))*y(4)^2 + 5*cos(y(3))^2*sin(y(1))*sin(y(5))*y(6)^2 + 5*sin(y(1))*sin(y(3))^2*sin(y(5))*y(6)^2))/(5*(sin(y(1))^2*sin(y(3))^2*sin(y(5))^2 + cos(y(1))^2*cos(y(3))^2*cos(y(5))^2 + 2*cos(y(1))^2*cos(y(3))^2*sin(y(5))^2 + 3*cos(y(1))^2*cos(y(5))^2*sin(y(3))^2 + 4*cos(y(3))^2*cos(y(5))^2*sin(y(1))^2 + 4*cos(y(1))^2*sin(y(3))^2*sin(y(5))^2 + 3*cos(y(3))^2*sin(y(1))^2*sin(y(5))^2 + 2*cos(y(5))^2*sin(y(1))^2*sin(y(3))^2 - 4*cos(y(1))*cos(y(3))*cos(y(5))^2*sin(y(1))*sin(y(3)) - 2*cos(y(1))*cos(y(3))^2*cos(y(5))*sin(y(1))*sin(y(5)) - 4*cos(y(1))*cos(y(3))*sin(y(1))*sin(y(3))*sin(y(5))^2 - 2*cos(y(1))*cos(y(5))*sin(y(1))*sin(y(3))^2*sin(y(5))))];
end