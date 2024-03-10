clear;
close all;

% midterm 23 Q5

Ds = tf([-1 -1],[3 2]);

Dz = c2d(Ds, 0.2, 'tustin');  