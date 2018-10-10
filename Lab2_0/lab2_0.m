%% LAB 2.0: SIGNALS AND SYSTEMS COMPUTATIONS USING MATLAB %%

%The goal of this lab is to gain familiarity with computations and plots 
%with Matlab, and to reinforce key concepts in signals and systems.  The 
%questions are chosen to illustrate how we can emulate continuous time 
%operations using the discrete time framework provided by Matlab.

clear;
clc;
close all;

%% FUNCTIONS AND PLOTS %%

t = -6:0.01:6;                          %Defining our time slot
x_b = signalx(t);                       %Calling the funtion for part b
x_c = signalx(t-3);                     %Calling the funtion for part c
x_d = signalx(3-t);                     %Calling the funtion for part d
x_e = signalx(2*t);                     %Calling the funtion for part e

figure(1);                              %Starting to plot
subplot(2,2,1);                         %First plot
plot(t, x_b, 'b');                      %Part b
title('Part b: signalx(t)');            %Plot title
xlabel('t');                            %X-axis is the time
subplot(2,2,2);                         %Second plot
plot(t, x_c, 'c');                      %Part c
title('Part c: signalx(t-3)');          %Plot title
xlabel('t');                            %X-axis is the time
subplot(2,2,3);                         %Third plot
plot(t, x_d, 'g');                      %Part d
title('Part d: signalx(3-t)');          %Plot title
xlabel('t');                            %X-axis is the time
subplot(2,2,4);                         %Fourth plot
plot(t, x_e, 'y');                      %Part e
title('Part e: signalx(2t)');           %Plot title
xlabel('t');                            %X-axis is the time
print('exercice_1','-dpng');            %Saving the plot

%% CONVOLUTION %%

dt = 0.01;                              %Sample spacing
t1 = -4 : dt : 4;
s1 = -2 : dt : -1;                      %Sampling times over [-2,-1]
s2 = 1 : dt : 3;                        %Sampling times over [1,3]
x1 = 3 * ((t1>=-2).*(t1<=-1));          %Given signal
x2 = 4 * ((t1>=1).*(t1<=3));            %Given signal


[y,t]= contconv(x1,x2,t1(1),t1(1),dt);  %My function
h = conv(x1, x2);                       %Built-in function
h = h/100;                              %Scaling

figure(2);                              %Starting to plot
subplot(2,1,1);                         %First plot
plot(t, y, 'b'); axis tight;            %Exercice function
title('Exercice function');             %Title
subplot(2,1,2);                         %Second plot
plot(h, 'r'); axis tight;               %Built-in function
title('MatLab Built-in function');      %Title
print('exercice_2','-dpng');            %Saving the plot

%% MATCH FILTER %%

dt = 0.01;                              %Sample spacing
t = -5 : dt : 5;                        %Time vector
s1 = 1 : dt : 3;                        %Sampling times over [1,3]
s2 = 2 : dt : 4;                        %Sampling times over [2,4]

u = 2 * ((t>=1).*(t<=3)) - 3 * ((t>=2).*(t<=4));	%Given signal
w = conj(fliplr(u));                                %Match filter

figure(3);                              %Starting to plot
subplot(2,1,1);                         %First plot
plot(t, u, 'r'); axis tight;            %Given signal
title('u(t))');                         %Title
subplot(2,1,2);                         %Second plot
plot(t, w, 'b'); axis tight;            %Match filter
title('u*(-t)');                        %Title
print('exercice_3a','-dpng');           %Saving the plot

[y,x]= contconv(u,w,t(1),t(1),dt);      %Convolution

figure(4);                              %Starting to plot
plot(x, y, 'b'); axis tight;            %Convolution
title('Convolution u(t) and u*(-t)');   %Title
print('exercice_3b','-dpng');           %Saving the plot

v = 2 * ((t>=-1).*(t<=2)) - 3 * ((t>=0).*(t<=1));	%Given signal
s = u + 1i*v;                                       %Total signal
d = -conj(s);                                       %Match filter

s_re = real(s);                         %Real part of the signal
d_re = real(d);                         %Real part of the match filter
s_im = imag(s);                         %Imaginary part of the signal
d_im = imag(d);                         %Imaginary part of the match filter

figure(5);                              %Starting to plot
subplot(2,2,1);                         %First plot
plot(t, s_re, 'r'); axis tight;         %Re[s(t)]
title('Real Part of s(t)');             %Title           
subplot(2,2,2);                         %Second plot
plot(t, d_re, 'b'); axis tight;         %Re[s*(-t)]
title('Real Part of s*(-t)');           %Title           
subplot(2,2,3);                         %Third plot
plot(t, s_im, 'r'); axis tight;         %Im[s(t)]
title('Imaginary Part of s(t)');        %Title 
subplot(2,2,4);                         %Fourth plot
plot(t, d_im, 'b'); axis tight;         %Im[s*(-t)]
title('Imaginary Part of s*(-t)');      %Title 
print('exercice_3c','-dpng');           %Saving the plot

[r,f]= contconv(s,d,t(1),t(1),dt);      %Convolution

r_re = real(r);                         %Real part of the convolution
r_im = imag(r);                         %Imaginary part of the convolution
r_ma = abs(r);                          %Magnitude of the convolution

figure(6);                                  %Starting to plot
subplot(3,1,1);                             %First plot
plot(f, r_re, 'r'); axis tight;             %Re[s(t)]
title('Real Part of the convolution');      %Title           
subplot(3,1,2);                             %First plot
plot(f, r_im, 'b'); axis tight;             %Re[s*(-t)]
title('Imaginary Part of the convolution'); %Title         
subplot(3,1,3);                             %First plot
plot(f, r_ma, 'c'); axis tight;             %Im[s(t)]
title('Magnitude of the convolution');      %Title  
print('exercice_3d','-dpng');               %Saving the plot

u1 = 2 * ((t>=-1).*(t<=1)) - 3 * ((t>=0).*(t<=2));      %Given signal
v1 = 2 * ((t>=-3).*(t<=0)) - 3 * ((t>=-2).*(t<=-1));	%Given signal
s1 = (u1 + 1i*v1) * exp(1j*pi/4);                       %Total signal

[g,h]= contconv(s1,d,t(1),t(1),dt);     %Convolution

g_re = real(g);                         %Real part of the convolution
g_im = imag(g);                         %Imaginary part of the convolution
g_ma = abs(g);                          %Magnitude of the convolution

figure(7);                                  %Starting to plot
subplot(3,1,1);                             %First plot
plot(h, g_re, 'r'); axis tight;             %Re[s(t)]
title('Real Part of the convolution');      %Title           
subplot(3,1,2);                             %First plot
plot(h, g_im, 'b'); axis tight;             %Re[s*(-t)]
title('Imaginary Part of the convolution'); %Title         
subplot(3,1,3);                             %First plot
plot(h, g_ma, 'c'); axis tight;             %Im[s(t)]
title('Magnitude of the convolution');      %Title  
print('exercice_3e','-dpng');               %Saving the plot

%% FOURIER TRANSFORM %%

dt = 1/16000;                           %Sample spacing
df_i = 1;                               %Desired frequency resolution
t = -8 : dt : 8;                        %Time vector
s = 3 * sinc(2*t - 3);                  %Input signal
[X, f, df] = contFT(s,t(1),dt,df_i);    %Given funtion
X_ma = abs(X);                          %Magnitude of the FT

figure(8);                              %Starting to plot
plot(f, X_ma, 'r'); axis tight;         %Starting to plot
title('Magnitude of the FT');           %Title           
print('exercice_4a','-dpng');           %Saving the plot

figure(9);                                      %Starting to plot
plot(f, X_ma, 'r'); xlim([-5,5]);               %Starting to plot
title('Close-up of the magnitud of the FT');    %Title           
print('exercice_4a2','-dpng');                  %Saving the plot

figure(10);                                     %Starting to plot
plot(f, X_ma, 'r'); xlim([-1.1,1.1]);           %Starting to plot
title('Close-up of the magnitud of the FT');    %Title           
print('exercice_4a3','-dpng');                  %Saving the plot

