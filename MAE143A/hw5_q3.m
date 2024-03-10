clear;
close all;

% Q3 
omega = pi;
t = 0:0.001:5;

yt = sin(omega*t);
%yt1 = tf(omega,[1 0 omega^2]);


h = 0.1;
t_k = 0:h:5.1; 
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

%y_foh = c2d(yt1,0.1,'foh');

figure(1);
hold on
plot(t,yt,'--b');
plot(t_k,y_k,'xk');
plot(t,y_zoh(1:5001),'-k');
%step(y_zoh)
xlim([0 5])
ylim([-1 1])

% h = 0.1

y_foh1 = t.*0;
counter = 2;

for i = 2:50 
    for j = 0.001:0.001:0.1
        mult = j == 0.10;
        mult2 = 0 == mult;
        if mult2 == 1
            y_foh1(counter)=(y_k(i)+(0.1)*(y_k(i)-y_k(i-1))/h);
        end
        if mult == 1
            y_foh1(counter) = y_k(i);
        end
        counter = counter+1;
    end
end



% h=0.2
% h2 = 0.2;
% t_k2 = 0:h2:5; 
% y_k2 = sin(omega*t_k2);
% 
% y_zoh2 = t.*0;
% counter = 1;

% for i = 1:25   
%     for j = 0.001:0.001:0.1
%     mult = j == 0.2;
%     mult2 = 0 == mult;
%     y_zoh(counter) = y_k(i)*mult2+y_k(i+1)*mult;
%     counter = counter+1;
%     end
% end
% 
% %h=0.4
% h4 = 0.4;
% t_k4 = 0:h4:5; 
% y_k4 = sin(omega*t_k4);
% 
% y_zoh4 = t.*0;
% counter = 1;
% 
% for i = 1:13   
%     for j = 0.001:0.001:0.1
%     mult = j == 0.4;
%     mult2 = 0 == mult;
%     y_zoh(counter) = y_k(i)*mult2+y_k(i+1)*mult;
%     counter = counter+1;
%     end
% end

figure(2);
hold on
plot(t,yt,'--b');
plot(t_k,y_k,'xk');
%plot(t_k2,y_k2,'xg');
%plot(t_k4,y_k4,'xc');
plot(t,y_foh1(1:5001),'-g')