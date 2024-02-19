clear;
close all;



%% 1

% A


figure(1);
FAP1 = RR_tf([1 -1],[1 1]); 
RR_bode(FAP1);


%% 2
%syms w
w = 0.1;

polycoeff1 = [1 0.77*w w];
polycoeff2 = [1 1.85*w w];

poles1 = RR_roots(polycoeff1);
poles2 = RR_roots(polycoeff2);

figure(2);
F2a = RR_tf( w^2,[1 0.77*w w^2]); 
RR_bode(F2a);

figure(3);
F2b = RR_tf( w^2,[1 1.85*w w^2]); 
RR_bode(F2b);

figure(4);
sys = tf(w^4,[1 (2.62*w) (2.4245*w^2+w^2) 2.62*w^3 w^4]);
step(sys)
