clear;
close all;

%% 1 2012
syms s M K C m k c

L1 = s^2*M + C*s + K;
L2 = C*s + K;
L3 = s^2*m + c*s + k;
L4 = c*s + k;
L5 = c*s + k;
L6 = L3*L1 - L5*L4 - L3*L5;
L7 = L3 * L2;

%% 2
zeta=0.1;
figure(1);
G2 = RR_tf([1 1 100],[100 200*zeta+10 20*zeta+100 10]); 
RR_bode(G2);



%% 3

figure(2);
G3 = RR_tf([100 -100],[1 0 100 0]);
RR_bode(G3);


%% 4
K = 4000*pi^2;
figure(3);
Y3 = RR_tf(K , [1000 0 -K]); % for W(s)=1
RR_bode(Y3);

figure(4);
Y32 = RR_tf(K , [1000 0 -K 0]); % for W(s)=1/s
RR_bode(Y32);


numcoeff1 = [1000 0 -K];
numcoeff2 = [1000 0 -K 0];
poles1 = RR_roots(numcoeff1);
poles2 = RR_roots(numcoeff2);
 
figure(5);
Y33 = tf(K, [1000 0 K]);
step(Y33);
axis([0 10 -10 10])


A1 = K/(poles1(1) - poles1(2));
A2 = K/((poles2(3))*(poles2(3) - poles2(1)));

B1 = K/(poles1(2) - poles1(1));
B2 = K/((poles2(1))*(poles2(1) - poles2(3)));

C2 = K/((poles2(2) - poles2(3))*(poles2(2) - poles2(1)));

