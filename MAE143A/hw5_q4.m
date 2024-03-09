clear;
close all;


Gs = RR_tf([-3 -2] ,[1 5 4 0]);

polycoeff = [1 5 4 0];

poles1 = RR_roots(polycoeff);

A = (-3*poles1(2)-2)/((poles1(2)+4)*(poles1(2)));
B = (-3*poles1(3)-2)/((poles1(3)+1)*(poles1(3)));
C = (-3*poles1(1)-2)/((poles1(1)+1)*(poles1(1)+4));