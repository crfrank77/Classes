clear;
close all;


% 2014 midterm 1
% q2

% A
F2=tf([1 0] ,[1 1]); 

h = 0.2;
[F2z] = c2d(F2,h,'tustin');

% B 
F1 = tf(1,[1 2 2 1]);  

F = F1*F2;

% C
[Faz] = c2d(F1,h,'tustin');

Fz = Faz*F2z;

F1b = RR_tf(1,[1 2 2 1]);  
F2b = RR_tf([1 0] ,[1 1]); 
F1bz = RR_C2D_tustin(F1b,h);
F2bz = RR_C2D_tustin(F2b,h);

Fzb = F1bz*F2bz;

[p,d,k,n] = RR_partial_fraction_expansion(Fzb);

