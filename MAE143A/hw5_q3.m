clear;
close all;

omega = pi;
t = 0:0.001:5;

yt = sin(omega*t);

h = 0.1;
t_k = 0:h:5; 
y_k = sin(omega*t_k);

y_zoh = t.*0;
counter = 1;

for i = 1:50   
    for j = 0.001:0.001:0.1
    mult = j == 0.1;
    mult2 = 0 == mult;
    y_zoh(counter) = y_k(i)*mult2+y_k(i+1)*mult;
    counter = counter+1;
    end
end

figure(1);
hold on
plot(t,yt,'--b');
plot(t_k,y_k,'xk');
plot(t,y_zoh(1:5001),'-k');