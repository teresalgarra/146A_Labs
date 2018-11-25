%% LAB 6.1: LINEAR MODULATION WITH TWO-DIMENSIONAL CONSTELLATIONS %%

%This is a follow-on to Software Lab 4.1, the code from which is our 
%starting point here. The objective is to implement in complex baseband a 
%linearly modulated system for a variety of signal constellations. We wish 
%to estimate the performance of these schemes for an ideal channel via 
%simulation, and to compare with analytical expressions. As in Software 
%Lab4.1, we use a trivial channel filter in this lab.

clear;
clc;
close all;

%% LABORATORY ASSIGNMENT %%

clear;
clc;
close all;

%% BPSK %%
a = 0.22;                                   %Bandwidth Excess
m = 4;                                      %Oversampling factor
l = 5;                                      %Length
[g, ~] = srrc_pulse(a, m, l);               %Using the created function

ns = 12000;                                 %Number of samples
x1 = randbit(ns);                           %Generating Random Bits
y = bpskmap(x1);                            %Mapping to BSPK

y1_u = 1+(ns-1)*m;                          %Upsampling by m
s_u = zeros(y1_u,1);                        %Upsampling by m                   
s_u(1:m:y1_u)=y*m;                          %Upsampling by m
tx_o = conv(s_u,g)/m;                       %Noiseless modulated signal
l_t = length(tx_o);                         %Length of tx_o
%t1 = cumsum(ones(l_t,1)/m)-1/m-5;          %Time vector for tx_o

SNR = 2.71;                                 %SNR from the Q function
Es = sum(abs(g).^2);                        %Energy per symbol
Eb = Es/1;                                  %Energy per bit
N0 = Eb/SNR;                                %Noice variance
noise = randn(l_t,1) + 1j.*randn(l_t,1);    %Noise vector
tx_o_n = tx_o + sqrt(N0/2).*noise;          %Noise added

% figure;                                     %Plotting
% plot(t1, tx_o_n, 'b'); axis tight;          %Function
% hold on;                                    %Holding on
% plot(t1, tx_o, 'r'); axis tight;            %Function
% title('Comparison with and without noise'); %Title
% xlabel('Samples');                          %Title for the x-axis
% print('images/exercise5_2','-dpng');        %Saving the plot

h = contconv(tx_o_n, g, 0, -5, 1/m);        %Sending b[n] through the system
t2 = cumsum(ones(length(h),1)/m)-1/m-10;    %Time vector for h

t0_index = find(t2==0);                     %Looking for the first zero
index0 = t0_index:m:m*ns+t0_index-1;        %Vector and steps
r = h(index0);                              %Getting the integers only
x2 = bpskmapinverse(r);                     %Reversing the mapping

error1 = 100 - sum(x2.'==x1)/ns*100;        %Error of the recovery

figure;                                     %Plotting
subplot(2,1,1);                             %First plot
plot(real(y),imag(y), 'bo'); axis tight;    %Function
title('Noiseless Constellation');           %Title
subplot(2,1,2);                             %First plot
plot(real(r),imag(r), 'co'); axis tight;    %Function
title('Constellation with Noise');          %Title
print('images/BPSK','-dpng');               %Saving the plot

%% 4PAM %%

a = 0.22;                                   %Bandwidth Excess
m = 4;                                      %Oversampling factor
l = 5;                                      %Length
[g, ~] = srrc_pulse(a, m, l);               %Using the created function

ns = 6000;                                  %Number of samples
x1 = randbit(ns);                           %Generating Random Bits
x2 = randbit(ns);                           %Generating Random Bits
y = fourpammap(x1,x2);                      %Mapping to 4PAM

y1_u = 1+(ns-1)*m;                          %Upsampling by m
s_u = zeros(y1_u,1);                        %Upsampling by m                   
s_u(1:m:y1_u)=y*m;                          %Upsampling by m
tx_o = conv(s_u,g)/m;                       %Noiseless modulated signal
l_t = length(tx_o);                         %Length of tx_o

SNR = 6.77;                                 %SNR from the Q function
Es = 5*sum(abs(g).^2);                      %Energy per symbol
Eb = Es/2;                                  %Energy per bit
N0 = Eb/SNR;                                %Noice variance
noise = randn(l_t,1) + 1j.*randn(l_t,1);    %Noise vector
tx_o_n = tx_o + sqrt(N0/2).*noise;          %Noise added

h = contconv(tx_o_n, g, 0, -5, 1/m);        %Sending b[n] through the system
t2 = cumsum(ones(length(h),1)/m)-1/m-10;    %Time vector for tx_o

t0_index = find(t2==0);                     %Looking for the first zero
index0 = t0_index:m:m*ns+t0_index-1;        %Vector and steps
r = h(index0);                              %Getting the integers only
[x12, x22] = fourpammapinverse(r);          %Reversing the mapping

error21 = 100 - sum(x12.'==x1)/ns*100;      %Error of the recovery
error22 = 100 - sum(x22.'==x2)/ns*100;      %Error of the recovery
error2 = (error21 + error22)./2;            %Error of recovery

figure;                                     %Plotting
subplot(2,1,1);                             %First plot
plot(real(y),imag(y), 'bo'); axis tight;    %Function
title('Noiseless Constellation');           %Title
subplot(2,1,2);                             %First plot
plot(real(r),imag(r), 'co'); axis tight;    %Function
title('Constellation with Noise');          %Title
print('images/4PAM','-dpng');               %Saving the plot

%% QPSK %%

a = 0.22;                                   %Bandwidth Excess
m = 4;                                      %Oversampling factor
l = 5;                                      %Length
[g, t] = srrc_pulse(a, m, l);               %Using the created function

ns = 6000;                                  %Number of samples
x1 = randbit(ns);                           %Generating Random Bits
x2 = randbit(ns);                           %Generating Random Bits
y = qpskmap(x1,x2);                         %Mapping to QPSK

y1_u = 1+(ns-1)*m;                          %Upsampling by m
s_u = zeros(y1_u,1);                        %Upsampling by m                   
s_u(1:m:y1_u)=y*m;                          %Upsampling by m
tx_o = conv(s_u,g)/m;                       %Noiseless modulated signal
l_t = length(tx_o);                         %Length of tx_o

SNR = 2.71;                                 %SNR from the Q function
Es = 2*sum(abs(g).^2);                      %Energy per symbol
Eb = Es/2;                                  %Energy per bit
N0 = Eb/SNR;                                %Noice variance
noise = randn(l_t,1) + 1j.*randn(l_t,1);    %Noise vector
tx_o_n = tx_o + sqrt(N0/2).*noise;          %Noise added

h = contconv(tx_o_n, g, 0, -5, 1/m);        %Sending b[n] through the system
t2 = cumsum(ones(length(h),1)/m)-1/m-10;    %Time vector for tx_o

t0_index = find(t2==0);                     %Looking for the first zero
index0 = t0_index:m:m*ns+t0_index-1;        %Vector and steps
r = h(index0);                              %Getting the integers only
[x12, x22] = qpskmapinverse(r);             %Reversing the mapping

error31 = 100 - sum(x12.'==x1)/ns*100;      %Error of the recovery
error32 = 100 - sum(x22.'==x2)/ns*100;      %Error of the recovery
error3 = (error31 + error32)./2;            %Error of recovery

figure;                                     %Plotting
subplot(2,1,1);                             %First plot
plot(real(y),imag(y), 'bo'); axis tight;    %Function
title('Noiseless Constellation');           %Title
subplot(2,1,2);                             %First plot
plot(real(r),imag(r), 'co'); axis tight;    %Function
title('Constellation with Noise');          %Title
print('images/QPSK','-dpng');               %Saving the plot

%% 16QAM %%

a = 0.22;                                   %Bandwidth Excess
m = 4;                                      %Oversampling factor
l = 5;                                      %Length
[g, ~] = srrc_pulse(a, m, l);               %Using the created function

ns = 3000;                                  %Number of samples
x1 = randbit(ns);                           %Generating Random Bits
x2 = randbit(ns);                           %Generating Random Bits
x3 = randbit(ns);                           %Generating Random Bits
x4 = randbit(ns);                           %Generating Random Bits
y = sixteenqammap(x1,x2,x3,x4);             %Mapping to 16QAM

y1_u = 1+(ns-1)*m;                          %Upsampling by m
s_u = zeros(y1_u,1);                        %Upsampling by m                   
s_u(1:m:y1_u)=y*m;                          %Upsampling by m
tx_o = conv(s_u,g)/m;                       %Noiseless modulated signal
l_t = length(tx_o);                         %Length of tx_o

SNR = 6.77;                                 %SNR from the Q function
Es = 10*sum(abs(g).^2);                     %Energy per symbol
Eb = Es/4;                                  %Energy per bit
N0 = Eb/SNR;                                %Noice variance
noise = randn(l_t,1) + 1j.*randn(l_t,1);    %Noise vector
tx_o_n = tx_o + sqrt(N0/2).*noise;          %Noise added

h = contconv(tx_o_n, g, 0, -5, 1/m);        %Sending b[n] through the system
t2 = cumsum(ones(length(h),1)/m)-1/m-10;    %Time vector for tx_o

t0_index = find(t2==0);                     %Looking for the first zero
index0 = t0_index:m:m*ns+t0_index-1;        %Vector and steps
r = h(index0);                              %Getting the integers only
[x12,x22,x32,x42]=sixteenqammapinverse(r);  %Reversing the mapping

error41 = 100 - sum(x12.'==x1)/ns*100;              %Error of the recovery
error42 = 100 - sum(x22.'==x2)/ns*100;              %Error of the recovery
error43 = 100 - sum(x32.'==x3)/ns*100;              %Error of the recovery
error44 = 100 - sum(x42.'==x4)/ns*100;              %Error of the recovery
error4 = (error41 + error42 + error43 + error44)./4;%Error of recovery

figure;                                     %Plotting
subplot(2,1,1);                             %First plot
plot(real(y),imag(y), 'bo'); axis tight;    %Function
title('Noiseless Constellation');           %Title
subplot(2,1,2);                             %First plot
plot(real(r),imag(r), 'co'); axis tight;    %Function
title('Constellation with Noise');          %Title
print('images/16QAM1','-dpng');             %Saving the plot

%% 8PSK %%

a = 0.22;                                   %Bandwidth Excess
m = 4;                                      %Oversampling factor
l = 5;                                      %Length
[g, ~] = srrc_pulse(a, m, l);               %Using the created function

ns = 4000;                                  %Number of samples
x1 = randbit(ns);                           %Generating Random Bits
x2 = randbit(ns);                           %Generating Random Bits
x3 = randbit(ns);                           %Generating Random Bits
y = eightpskmap(x1,x2,x3);                  %Mapping to 8PSK

y1_u = 1+(ns-1)*m;                          %Upsampling by m
s_u = zeros(y1_u,1);                        %Upsampling by m                   
s_u(1:m:y1_u)=y*m;                          %Upsampling by m
tx_o = conv(s_u,g)/m;                       %Noiseless modulated signal
l_t = length(tx_o);                         %Length of tx_o

SNR = 6.16;                                 %SNR from the Q function
Es = sum(abs(g).^2);                        %Energy per symbol
Eb = Es/3;                                  %Energy per bit
N0 = Eb/SNR;                                %Noice variance
noise = randn(l_t,1) + 1j.*randn(l_t,1);    %Noise vector
tx_o_n = tx_o + sqrt(N0/2).*noise;          %Noise added

h = contconv(tx_o_n, g, 0, -5, 1/m);        %Sending b[n] through the system
t2 = cumsum(ones(length(h),1)/m)-1/m-10;    %Time vector for tx_o

t0_index = find(t2==0);                     %Looking for the first zero
index0 = t0_index:m:m*ns+t0_index-1;        %Vector and steps
r = h(index0);                              %Getting the integers only
[x12,x22,x32]=eightpskmapinverse(r);        %Reversing the mapping

error51 = 100 - sum(x12.'==x1)/ns*100;      %Error of the recovery
error52 = 100 - sum(x22.'==x2)/ns*100;      %Error of the recovery
error53 = 100 - sum(x32.'==x3)/ns*100;      %Error of the recovery
error5 = (error51 + error52 + error53)./3;  %Error of recovery

figure;                                     %Plotting
subplot(2,1,1);                             %First plot
plot(real(y),imag(y), 'bo'); axis tight;    %Function
title('Noiseless Constellation');           %Title
subplot(2,1,2);                             %First plot
plot(real(r),imag(r), 'co'); axis tight;    %Function
title('Constellation with Noise');          %Title
print('images/8PSK','-dpng');               %Saving the plot

%% 13A.- 16QAM %%

ns = 3000;                                  %Number of samples
x1 = randbit(ns);                           %Generating Random Bits
x2 = randbit(ns);                           %Generating Random Bits
x3 = randbit(ns);                           %Generating Random Bits
x4 = randbit(ns);                           %Generating Random Bits
y = sixteenqammap(x1,x2,x3,x4);             %Mapping to 16QAM

SNR = 6.77;                                 %SNR from the Q function
Es = 10;                                    %Energy per symbol
Eb = Es/4;                                  %Energy per bit
N0 = Eb/SNR;                                %Noice variance
noise = randn(ns,1) + 1j.*randn(ns,1);      %Noise vector
r = y + sqrt(N0/2).*noise.';                %Noise added

[x12,x22,x32,x42]=sixteenqammapinverse(r);  %Reversing the mapping

error61 = 100 - sum(x12.'==x1)/ns*100;              %Error of the recovery
error62 = 100 - sum(x22.'==x2)/ns*100;              %Error of the recovery
error63 = 100 - sum(x32.'==x3)/ns*100;              %Error of the recovery
error64 = 100 - sum(x42.'==x4)/ns*100;              %Error of the recovery
error6 = (error61 + error62 + error63 + error64)./4;%Error of recovery

figure;                                     %Plotting
subplot(2,1,1);                             %First plot
plot(real(y),imag(y), 'bo'); axis tight;    %Function
title('Noiseless Constellation');           %Title
subplot(2,1,2);                             %First plot
plot(real(r),imag(r), 'co'); axis tight;    %Function
title('Constellation with Noise');          %Title
print('images/16QAM2','-dpng');             %Saving the plot

%% 13B.- 16QAM %%

ns = 3000;                                  %Number of samples
x1 = randbit(ns);                           %Generating Random Bits
x2 = randbit(ns);                           %Generating Random Bits
x3 = randbit(ns);                           %Generating Random Bits
x4 = randbit(ns);                           %Generating Random Bits
y = sixteenqammap(x1,x2,x3,x4);             %Mapping to 16QAM

SNR = 9.77;                                 %SNR from the Q function
Es = 10;                                    %Energy per symbol
Eb = Es/4;                                  %Energy per bit
N0 = Eb/SNR;                                %Noice variance
noise = randn(ns,1) + 1j.*randn(ns,1);      %Noise vector
r = y + sqrt(N0/2).*noise.';                %Noise added

[x12,x22,x32,x42]=sixteenqammapinverse(r);  %Reversing the mapping

error71 = 100 - sum(x12.'==x1)/ns*100;              %Error of the recovery
error72 = 100 - sum(x22.'==x2)/ns*100;              %Error of the recovery
error73 = 100 - sum(x32.'==x3)/ns*100;              %Error of the recovery
error74 = 100 - sum(x42.'==x4)/ns*100;              %Error of the recovery
error7 = (error71 + error72 + error73 + error74)./4;%Error of recovery

figure;                                     %Plotting
subplot(2,1,1);                             %First plot
plot(real(y),imag(y), 'bo'); axis tight;    %Function
title('Noiseless Constellation');           %Title
subplot(2,1,2);                             %First plot
plot(real(r),imag(r), 'co'); axis tight;    %Function
title('Constellation with Noise');          %Title
print('images/16QAM3','-dpng');             %Saving the plot
