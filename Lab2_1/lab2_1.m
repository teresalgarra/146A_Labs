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
dt = 0.01;
osfactor = 1/dt;
N = 100;
tt = dt:dt:N;

%%1.1: u_c and u_s%%

u_c = [];
for n = 1:1:N
    pos_c = randi(length(set));           %This selects a random position of the set for b_c
    b_c = set(pos_c);                     %We get the coefficient b_c
    bcc = repmat(b_c, osfactor, 1);
    bbc = bcc(:)';
    u_c = [u_c, bbc];
end

u_s = [];
for n = 1:1:N
    pos_s = randi(length(set));           %This selects a random position of the set for b_c
    b_s = set(pos_s);                     %We get the coefficient b_c
    bss = repmat(b_s, osfactor, 1);
    bbs = bss(:)';
    u_s = [u_s, bbs];
end

figure;
subplot(2,1,1);
plot(tt, u_c, 'b'); xlim([1,10]);
title('Uc');
xlabel('t');
subplot(2,1,2);
plot(tt, u_s, 'r'); xlim([1,10]);
title('Us');
xlabel('t');
print('images/1_1','-dpng');

%%1.2: Upconversion%%


u_p1 = u_c.*cos(40.*pi.*tt);            %We apply the formula

figure;
plot(tt, u_p1, 'b'); xlim([1,4]);
title('Up without Q');
xlabel('t');
print('images/1_2','-dpng');

%%1.3: Passband signal%%

u_p = u_p1 - u_s.*sin(40*pi*tt);

figure;
plot(tt, u_p, 'r'); xlim([1,4]);
title('Up with Q');
xlabel('t');
print('images/1_3','-dpng');

%%1.4: Downconversion%%

lpf = ((tt>=0).*(tt<=0.25));

v_c1 = u_p.*2.*cos(40*pi*tt);
[v_c, t_c] = contconv(v_c1,lpf,0,0,0.01);

v_s1 = u_p.*2.*sin(40*pi*tt);
[v_s, t_s] = contconv(v_s1,lpf,0,0,0.01);

figure;
subplot(2,1,1);
plot(t_c, v_c, 'b'); xlim([1,N]);
title('Vc');
xlabel('t');
subplot(2,1,2);
plot(t_s, v_s, 'r'); xlim([1,N]);
title('Vs');
xlabel('t');
print('images/1_4','-dpng');

%%1.5: Downconversion with phase%%

w_c1 = u_p.*2.*cos(40*pi*tt + pi/4);
[w_c, tt_c] = contconv(w_c1,lpf,0,0,0.01);

w_s1 = u_p.*2.*sin(40*pi*tt + pi/4);
[w_s, tt_s] = contconv(w_s1,lpf,0,0,0.01);

figure;
subplot(2,1,1);
plot(tt_c, w_c, 'b'); xlim([1,N]);
title('Vc');
xlabel('t');
subplot(2,1,2);
plot(tt_s, w_s, 'r'); xlim([1,N]);
title('Vs');
xlabel('t');
print('images/1_5','-dpng');

%%1.6: Recovering uc and us%%

u_c2 = real(w_c * exp(1j*pi/4));
u_s2 = real(w_s * exp(1j*pi/4));

figure;
subplot(2,1,1);
plot(tt_c, u_c2, 'b'); xlim([1,N]);
title('Recovered Uc');
xlabel('t');
subplot(2,1,2);
plot(tt_s, u_s2, 'r'); xlim([1,N]);
title('Recovered Us');
xlabel('t');
print('images/1_6','-dpng');

figure;
subplot(2,1,1);
plot(tt_c, u_c2, 'b'); xlim([1,N]);
title('Recovered Uc');
xlabel('t');
subplot(2,1,2);
plot(tt, u_c, 'r'); xlim([1,N]);
title('Uc');
xlabel('t');
print('images/1_6a','-dpng');

figure;
subplot(2,1,1);
plot(tt_c, u_s2, 'b'); xlim([1,N]);
title('Recovered Us');
xlabel('t');
subplot(2,1,2);
plot(tt, u_s, 'r'); xlim([1,N]);
title('Us');
xlabel('t');
print('images/1_6b','-dpng');

