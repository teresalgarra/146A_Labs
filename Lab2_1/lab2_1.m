%% LAB 2.1: MODELING CARRIER PHASE UNCERTAINTY %%

%The goal of this lab is to explore modeling and receiver operations in 
%complex baseband. In particular, we model and undo the effect of carrier 
%phase mismatch between the receiver LO and the incoming carrier.

clear;
clc;
close all;

%% LABORATORY ASSIGNMENT %%

%%Random coeficients%%
set = [-1,1];                         %These are the possible coefficients
pos_c = randi(length(set));           %This selects a random position of the set for b_c
b_c = set(pos_c);                     %We get the coefficient b_c
pos_s = randi(length(set));           %This selects a random position of the set for b_s
b_s = set(pos_s);                     %We get the coefficient b_s

%%Rectangular pulse%%
dt = 0.001;
N = 100;
tt = 1:dt:N;
tt = tt';
p = (tt>=0).*(tt<=1);                   %Given signal

l = length(tt);

%%1.1: u_c and u_s%%

symbols = 10;
u_c = zeros(l, 1);
for n = 1:1:l
    for t = 1:1:l
        for i = 1:1:symbols
            pos_c = randi(length(set));           %This selects a random position of the set for b_c
            b_c = set(pos_c);                     %We get the coefficient b_c
        end
        if t>n
            u_c(t) = b_c.*p(t-n);                 %We create the function u_c
        end
    end
end

u_s = zeros(l, 1);
for n = 1:1:l
    for t = 1:1:l
        for i = 1:1:symbols
            pos_s = randi(length(set));           %This selects a random position of the set for b_c
            b_s = set(pos_s);                     %We get the coefficient b_s
        end
        if t>n
            u_s(t) = b_s.*p(t-n);                 %We create the function u_s
        end
    end
end

figure;
subplot(2,1,1);
plot(tt, u_c, 'b'); axis tight;
title('Uc');
xlabel('t');
subplot(2,1,2);
plot(tt, u_s, 'r'); axis tight;
title('Us');
xlabel('t');
print('images/1_1','-dpng');

%%1.2: Upconversion%%

symbols = 4;
u_p1 = zeros(l, 1);
u_c = zeros(l, 1);
for n = 1:1:l
    for t = 1:1:l
        for i = 1:1:symbols
            pos_c = randi(length(set));           %This selects a random position of the set for b_c
            b_c = set(pos_c);                     %We get the coefficient b_c
        end
        if t>n
           u_c(t) = b_c.*p(t-n);                 %We create the function u_c
           u_p1(t) = u_c(t).*cos(40*pi*t);       %We apply the formula
        end
    end
end

figure;
plot(tt, u_p1, 'b'); axis tight;
title('Up without Q');
xlabel('t');
print('images/1_2','-dpng');

%%1.3: Passband signal%%

u_p = zeros(l, 1);
u_c = zeros(l, 1);
u_s = zeros(l, 1);
for n = 1:1:l
    for t = 1:1:l
        for i = 1:1:symbols
            pos_c = randi(length(set));           %This selects a random position of the set for b_c
            b_c = set(pos_c);                     %We get the coefficient b_c
            pos_s = randi(length(set));           %This selects a random position of the set for b_c
            b_s = set(pos_s);                     %We get the coefficient b_s
        end
        if t>n
            u_c(t) = b_c.*p(t-n);                 %We create the function u_c
            u_s(t) = b_s.*p(t-n);                %We create the function u_s
            u_p(t) = u_c(t).*cos(40*pi*t) - u_s(t).*sin(40*pi*t);   %We apply the formula
        end
    end
end

figure;
plot(tt, u_p, 'r'); axis tight;
title('Up with Q');
xlabel('t');
print('images/1_3','-dpng');

%%1.4: Downconversion%%

lpf = ((tt>=0).*(tt<=0.25));

v_c1 = 2.* u_p.*cos(40*pi*tt);
[v_c, t_c] = contconv(v_c1,lpf,0,0,0.01);

v_s1 = 2.*u_p.*sin(40*pi*tt);
[v_s, t_s] = contconv(v_s1,lpf,0,0,0.01);

figure;
subplot(2,1,1);
plot(t_c, v_c, 'b'); axis tight;
title('Vc');
xlabel('t');
subplot(2,1,2);
plot(t_s, v_s, 'r'); axis tight;
title('Vs');
xlabel('t');
print('images/1_4','-dpng');

%%1.5: Downconversion with phase%%

w_c1 = 2.* u_p.*cos(40*pi*tt + pi/4);
[w_c, tt_c] = contconv(w_c1,lpf,0,0,0.01);

w_s1 = 2.*u_p.*sin(40*pi*tt + pi/4);
[w_s, tt_s] = contconv(w_s1,lpf,0,0,0.01);

figure;
subplot(2,1,1);
plot(tt_c, w_c, 'b'); axis tight;
title('Vc');
xlabel('t');
subplot(2,1,2);
plot(tt_s, w_s, 'r'); axis tight;
title('Vs');
xlabel('t');
print('images/1_5','-dpng');

%%1.6: Recovering uc and us%%

u_c2 = w_c * exp(1j*pi/4);
u_s2 = w_s * exp(1j*pi/4);

figure;
subplot(2,1,1);
plot(tt_c, u_c2, 'b'); axis tight;
title('Recovered Uc');
xlabel('t');
subplot(2,1,2);
plot(tt_s, u_s2, 'r'); axis tight;
title('Recovered Us');
xlabel('t');
print('images/1_6','-dpng');

figure;
subplot(2,1,1);
plot(tt_c, u_c2, 'b'); axis tight;
title('Recovered Uc');
xlabel('t');
subplot(2,1,2);
plot(tt, u_c, 'r'); axis tight;
title('Uc');
xlabel('t');
print('images/1_6a','-dpng');

figure;
subplot(2,1,1);
plot(tt_c, u_2, 'b'); axis tight;
title('Recovered Us');
xlabel('t');
subplot(2,1,2);
plot(tt, u_s, 'r'); axis tight;
title('Us');
xlabel('t');
print('images/1_6b','-dpng');

