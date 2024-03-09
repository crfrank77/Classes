clear;
close all;


% 2014 midterm 1

polycoeff1 = [1 2 2 1 0];

poles1 = RR_roots(polycoeff1);


A = 1 / ((poles1(2)-poles1(3))*(poles1(2)-poles1(4))*(poles1(2)));

B = 1 / ((poles1(3)-poles1(2))*(poles1(3)-poles1(4))*(poles1(3)));

C = 1 / ((poles1(4)-poles1(2))*(poles1(4)-poles1(3))*(poles1(4)));

D = 1 / ((poles1(1)-poles1(2))*(poles1(1)-poles1(3))*(poles1(1)-poles1(4)));

Y=RR_tf(1,[1 2 2 1]); 
figure(1);
RR_step(Y);



% D 
h = 0.2;
[Yz] = RR_C2D_tustin(Y,h);

% E
omegac = 8;
[Yz2] = RR_C2D_tustin(Y,h,omegac);

figure(2);
RR_step(Yz);