clear;

%% 1 2023
syms s

lim10 = limit((-3*s-2)/((s+4)*(s+1)),s,Inf);
lim1inf = limit((-3*s-2)/((s+4)*(s+1)),s,0);
% Y(s)
polycoeff1 = [1 5 4 0];
numcoeff1 = [-3 -2];

roots1 = RR_roots(numcoeff1);
poles1 = RR_roots(polycoeff1);

% poles1(3) = -4
A = (numcoeff1(1)*poles1(3)+numcoeff1(2))/((poles1(3))*(poles1(3)+1));

% poles1(2) = -1
B = (numcoeff1(1)*poles1(2)+numcoeff1(2))/((poles1(2))*(poles1(2)+4));

%poles1(1) = 0
C = (numcoeff1(1)*poles1(1)+numcoeff1(2))/((poles1(1)+4)*(poles1(1)+1));

% plot
figure(1);
sys = tf([-3 -2],[1 5 4 0]);
step(sys)

%% 2 2023

% T(s)
polycoeff2 = [1 2 2 0];

roots2 = RR_roots(numcoeff1);
poles2 = RR_roots(polycoeff2); 

lim20 = limit((-3*s-2)/(s^2+2*s+2),s,Inf);
lim2inf = limit((-3*s-2)/(s^2+2*s+2),s,0);

% poles2(3) = -1+i
A2 = (numcoeff1(1)*poles2(3)+numcoeff1(2))/((poles2(3))*(poles2(3)-poles2(2)));

% poles2(2) = -1-i
B2 = (numcoeff1(1)*poles2(2)+numcoeff1(2))/((poles2(2))*(poles2(2)-poles2(3)));

%poles2(1) = 0
C2 = (numcoeff1(1)*poles2(1)+numcoeff1(2))/((poles2(1)-poles2(3))*(poles2(1)-poles2(2)));

%plot
figure(2);
sys2 = tf([-3 -2],[1 2 2 0]);
step(sys2)


%% 3 2023

sys3 = tf(1,[1 5 0]);

polycoeff3 = [1 5 0];

poles3 = RR_roots(polycoeff3); 

lim30 = limit((1)/(s+5),s,Inf);
lim3inf = limit((1)/(s+5),s,0);

% poels3(1) = 0
B3 = 1/(poles3(1)+5);

% poels3(2) = -5
A3 = 1/poles3(2);

%plot
figure(3);
sys3 = tf(1,[1 5 0]);
step(sys3)

