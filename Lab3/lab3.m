%% LAB 4.1: LINEAR MODULATION OVER A NOISELESS IDEAL CHANNEL %%

%This is the first of a sequence of software labs which gradually develop a
%reasonably complete Matlab simulator for a linearly modulated system. The 
%follow-on labs are Software Lab 6.1 in Chapter 6, and Software Lab 8.1 in 
%Chapter 8.

clear;
clc;
close all;

%% LABORATORY ASSIGNMENT %%

clear;
clc;
close all;

%%Exercice 1%%
a = 0.22;                                   %Bandwidth Excess
m = 4;                                      %Oversampling factor
l = 5;                                      %Length
[g, t] = srrc_pulse(a, m, l);               %Using the created function
figure;                                     %Starting to plot
plot(t, g, 'r'); axis tight;                %Function
title('SRRC Pulse');                        %Title
xlabel('t');                                %Title for the x-axis
print('images/exercise1','-dpng');          %Saving the plot

%%Exercice 2%%
[G,f,~] = contFT(g,-5,1/m,1/64);            %DFT of g(t)
G_ma = abs(G);                              %Magnitude of G(f)
figure;                                     %Starting to plot
plot(f, G_ma, 'b'); axis tight;             %Function
title('DFT of the SRRC Pulse');             %Title
xlabel('f');                                %Title for the x-axis
print('images/exercise2','-dpng');          %Saving the plot

%%Exercice 3%%
ns = 100;                                   %Number of symbols
v = 1:1:100;
a = randi([0,1],100,1);                     %Random coefficients
b = zeros(100,1);                           %Preallocating b[n]
for i = 1:100                               %Loop size
    if a(i) == 0                            %Mapping 0s to 1s
        b(i) = 1;                           %Mapping 0s to 1s
    elseif a(i) == 1                        %Mapping 1s to -1s
        b(i) = -1;                          %Mapping 1s to -1s
    end                                     %End of conditions
end                                         %End of loop
figure;                                     %Starting to plot
plot(v, b, 'r'); axis tight;                %Function
title('b[n]');                              %Title
xlabel('t');                                %Title for the x-axis
print('images/exercise3','-dpng');          %Saving the plot


%%Exercice 4%%
ns_u = 1+(ns-1)*m;                          %Upsampling by m
s_u = zeros(ns_u,1);                        %Upsampling by m                   
s_u(1:m:ns_u)=b;                            %Upsampling by m
tx_o = conv(s_u,g);                         %Noiseless modulated signal
t1 = cumsum(ones(length(tx_o),1)/m)-1/m-5;  %Time vector for tx_o
t2 = cumsum(ones(length(s_u),1)/m)-1/m;     %Time vector for s_u
figure;                                     %Starting to plot
plot(t1, tx_o, 'r'); axis tight;            %Function
hold on;                                    %Holding on the plot
plot(t2, s_u, 'b'); axis tight;             %Function
title('Convolution of b[n] and g(t)');      %Title
xlabel('t');                                %Title for the x-axis
print('images/exercise4','-dpng');          %Saving the plot


%%Exercice 5%%
h = contconv(tx_o, g, 0, -5, 1/m);          %Sending b[n] through the system
t3 = cumsum(ones(length(h),1)/m)-1/m-10;    %Time vector for tx_o
t4 = cumsum(ones(length(s_u),1)/m)-1/m;     %Time vector for s_u
figure;                                     %Starting to plot
plot(t3, h, 'r'); axis tight;               %Function
hold on;                                    %Holding on the plot
plot(t4, s_u, 'b'); axis tight;             %Function
title('Convolution of tx_o(t) and g(t)');   %Title
xlabel('t');                                %Title for the x-axis
print('images/exercise5','-dpng');          %Saving the plot

%%Exercice 6%%
t0_index = find(t1==0);                     %Looking for the first zero
index0 = t0_index:m:m*ns+t0_index;          %Vector and steps
r = tx_o(index0);                         	%Getting the integers only
d = zeros(100,1);                           %Preallocating d[n]
for i = 1:100                               %Loop size
    if r(i) < 0                             %Mapping 0s to -1s
        d(i) = 1;                           %Mapping 0s to -1s
    elseif r(i) > 0                         %Mapping 1s to 1s
        d(i) = 0;                           %Mapping 1s to 1s
    end                                     %End of conditions
end                                         %End of loop
figure;                                     %Starting to plot
subplot(2,1,1);                             %First plot
plot(v, d, 'r'); axis tight;                %Function
title('Recovering a[n] from tx_o(t)');      %Title
xlabel('t');                                %Title for the x-axis
subplot(2,1,2);                             %Second plot
plot(v, a, 'b'); axis tight;                %Function
title('a[n]');                              %Title
xlabel('t');                                %Title for the x-axis
print('images/exercise6','-dpng');          %Saving the plot

error1 = 100 - sum(d==a)/ns*100;            %Error of the recovery

%%Exercice 7%%

t1_index = find(t3==0);                     %Looking for the first zero
index1 = t1_index:m:m*ns+t1_index;          %Vector and steps
y = h(index1);                              %Getting the integers only
f = zeros(100,1);                           %Preallocating d[n]
for i = 1:100                               %Loop size
    if y(i) < 0                             %Mapping 0s to -1s
        f(i) = 1;                           %Mapping 0s to -1s
    elseif y(i) > 0                         %Mapping 1s to 1s
        f(i) = 0;                           %Mapping 1s to 1s
    end                                     %End of conditions
end                                         %End of loop
figure;                                     %Starting to plot
subplot(2,1,1);                             %First plot
plot(v, f, 'r'); axis tight;                %Function
title('Recovering a[n] from tx_o(t)');      %Title
xlabel('t');                                %Title for the x-axis
subplot(2,1,2);                             %Second plot
plot(v, a, 'b'); axis tight;                %Function
title('a[n]');                              %Title
xlabel('t');                                %Title for the x-axis
print('images/exercise7','-dpng');          %Saving the plot

error2 = 100 - sum(f==a)/ns*100;            %Error of the recovery

%%Exercice 8%%
df = 1/40;                                  %Given phase offset
u = tx_o.*exp(-1j*(pi*2*df.*t1+pi/2));      %Getting the phase offset
w = conv(u, g)/m;                           %Convolution with the filter
t2_index = find(t3==0);                     %Looking for the first zero
index2 = t2_index:m:m*ns+t2_index-1;        %Vector and steps
w = w(index2);                              %Getting the integers only
figure;                                     %Starting to plot
plot(real(w), imag(w), 'ro'); axis tight;   %Function
title('Imaginary part vs Real part');       %Title
xlabel('Real');                             %Title for the x-axis
ylabel('Imaginary');                        %Title for the x-axis
print('images/exercise8','-dpng');          %Saving the plot

%%Exercice 9%%
op = zeros(1,ns);                           %Preallocating the output
op(1) = b(1);                               %Initialiting the vector
for i=2:ns                                  %Loop
    op(i) = w(i).*conj(w(i-1)).*op(i-1);    %Output
    op(i) = sign(real(op(i)));              %Sign and real value of the output
end                                         %End of the loop
figure;                                     %Starting to plot
subplot(2,1,1);                             %First plot
plot(v, op, 'r'); axis tight;               %Function
title('Recovering b[n] from  the output');  %Title
xlabel('t');                                %Title for the x-axis
subplot(2,1,2);                             %Second plot
plot(v, b, 'b'); axis tight;                %Function
title('b[n]');                              %Title
xlabel('t');                                %Title for the x-axis
print('images/exercise9','-dpng');          %Saving the plot

error3 = 100 - sum(op'==b)/ns*100;          %Error of the recovery

%SECOND SYSTEM TRY
% c = zeros(ns+1,1);                          %Starting c[n]
% c(1) = 1;                                   %First value
% for n = 1:ns                                %Loop for the function
%     if a(n) == 0                            %First condition
%         c(n+1) = c(n);                      %Positive value of c(n)
%     elseif a(n) == 1                        %Second condition
%         c(n+1) = -c(n);                     %Negative value of c(n)
%     end                                     %End of conditions
% end                                         %End of loop
% c = c(2:end);                               %Removing the first value of c
% ns_uc = 1+(ns-1)*m;                         %Upsampling by m
% s_uc = zeros(ns_uc,1);                      %Upsampling by m                   
% s_uc(1:m:ns_uc)=c;                          %Upsampling by m
% c1 = conv(s_uc,g);                          %Noiseless modulated signal
% c2 = conv(c1,g);                            %Noiseless modulated signal
% t5 = cumsum(ones(length(c2),1)/m)-1/m-10;   %Time vector for c2
% t6 = cumsum(ones(length(s_uc),1)/m)-1/m;    %Time vector for s_uc
% 
% figure;                                     %Starting to plot
% plot(t5, c2, 'r'); axis tight;              %Function
% hold on;                                    %Holding on the plot
% plot(t6, s_uc, 'b'); axis tight;            %Function
% title('c[n] convolved with g(t) twice');    %Title
% xlabel('t');                                %Title for the x-axis
% print('images/exercise9','-dpng');          %Saving the plot
