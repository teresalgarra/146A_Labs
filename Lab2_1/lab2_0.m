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
print('images/exercice_1','-dpng');     %Saving the plot

%% CONVOLUTION %%

dt = 0.01;                              %Sample spacing
t = -4 : dt : 4;                        %Time vector
x1 = 3*((t>=-2).*(t<=-1));              %Given signal
x2 = 4*((t>=1).*(t<=3));                %Given signal

[y,x]= contconv(x1,x2,t(1),t(1),dt);    %My function

figure(2);                              %Starting to plot
plot(x, y, 'b'); axis tight;            %Exercice function
title('ContConv function');             %Title
xlabel('t');                            %Title for the x-axis
print('images/exercice_2','-dpng');     %Saving the plot

%% MATCH FILTER %%

dt = 0.01;                              %Sample spacing
t = -5 : dt : 5;                        %Time vector

u = 2*((t>=1).*(t<=3))-3*((t>=2).*(t<=4));	%Given signal
w = conj(fliplr(u));                        %Match filter

figure(3);                              %Starting to plot
subplot(2,1,1);                         %First plot
plot(t, u, 'b'); axis tight;            %Given signal
title('u(t))');                         %Title
xlabel('t');                            %Title for the x-axis
subplot(2,1,2);                         %Second plot
plot(t, w, 'c'); axis tight;            %Match filter
title('u*(-t)');                        %Title
xlabel('t');                            %Title for the x-axis
print('images/exercice_3a','-dpng');    %Saving the plot

[y,x]= contconv(u,w,t(1),t(1),dt);      %Convolution

figure(4);                              %Starting to plot
plot(x, y, 'g'); axis tight;            %Convolution
title('Convolution u(t) and u*(-t)');   %Title
xlabel('t');                            %Title for the x-axis
print('images/exercice_3b','-dpng');    %Saving the plot

v = 2*((t>=-1).*(t<=2))-3*((t>=0).*(t<=1));	%Given signal
s = u + 1i*v;                               %Total signal
d = fliplr(conj(s));                        %Match filter

s_re = real(s);                         %Real part of the signal
d_re = real(d);                         %Real part of the match filter
s_im = imag(s);                         %Imaginary part of the signal
d_im = imag(d);                         %Imaginary part of the match filter

figure(5);                              %Starting to plot
subplot(2,2,1);                         %First plot
plot(t, s_re, 'b'); axis tight;         %Re[s(t)]
title('Real Part of s(t)');             %Title
xlabel('t');                            %Title for the x-axis
subplot(2,2,2);                         %Second plot
plot(t, d_re, 'c'); axis tight;         %Re[s*(-t)]
title('Real Part of s*(-t)');           %Title
xlabel('t');                            %Title for the x-axis
subplot(2,2,3);                         %Third plot
plot(t, s_im, 'g'); axis tight;         %Im[s(t)]
title('Imaginary Part of s(t)');        %Title
xlabel('t');                            %Title for the x-axis
subplot(2,2,4);                         %Fourth plot
plot(t, d_im, 'y'); axis tight;         %Im[s*(-t)]
title('Imaginary Part of s*(-t)');      %Title
xlabel('t');                            %Title for the x-axis
print('images/exercice_3c','-dpng');    %Saving the plot

[r,f]= contconv(s,d,t(1),t(1),dt);      %Convolution

r_re = real(r);                         %Real part of the convolution
r_im = imag(r);                         %Imaginary part of the convolution
r_ma = abs(r);                          %Magnitude of the convolution

figure(6);                                  %Starting to plot
subplot(3,1,1);                             %First plot
plot(f, r_re, 'b'); axis tight;             %Re[s(t)]
title('Real Part of the convolution');      %Title
xlabel('t');                                %Title for the x-axis
subplot(3,1,2);                             %Second plot
plot(f, r_im, 'c'); axis tight;             %Re[s*(-t)]
title('Imaginary Part of the convolution'); %Title
xlabel('t');                                %Title for the x-axis
subplot(3,1,3);                             %Third plot
plot(f, r_ma, 'g'); axis tight;             %Im[s(t)]
xlabel('t');                                %Title for the x-axis
title('Magnitude of the convolution');      %Title
print('images/exercice_3d','-dpng');        %Saving the plot

u1 = 2 * ((t>=-1).*(t<=1)) - 3 * ((t>=0).*(t<=2));      %Given signal
v1 = 2 * ((t>=-3).*(t<=0)) - 3 * ((t>=-2).*(t<=-1));	%Given signal
s1 = (u1 + 1i*v1) * exp(1j*pi/4);                       %Total signal

[g,h]= contconv(s1,d,t(1),t(1),dt);         %Convolution

g_re = real(g);                             %Real part of the convolution
g_im = imag(g);                             %Imaginary part of the convolution
g_ma = abs(g);                              %Magnitude of the convolution

figure(7);                                  %Starting to plot
subplot(3,1,1);                             %First plot
plot(h, g_re, 'b'); axis tight;             %Re[s(t)]
title('Real Part of the convolution');      %Title
xlabel('t');                                %Title for the x-axis
subplot(3,1,2);                             %Second plot
plot(h, g_im, 'c'); axis tight;             %Re[s*(-t)]
title('Imaginary Part of the convolution'); %Title
xlabel('t');                                %Title for the x-axis
subplot(3,1,3);                             %Third plot
plot(h, g_ma, 'g'); axis tight;             %Im[s(t)]
title('Magnitude of the convolution');      %Title
xlabel('t');                                %Title for the x-axis
print('images/exercice_3e','-dpng');        %Saving the plot

%% FOURIER TRANSFORM %%

dt = 1/(16e6);                          %Sample spacing
df_i = 1e3;                             %Desired frequency resolution
t = -8e-6 : dt : 8e-6;                  %Time vector
s = 3 * sinc(2*t*1e6 - 3);              %Input signal
[X, f, ~] = contFT(s,t(1),dt,df_i);     %Given funtion
X_ma = abs(X);                          %Magnitude of the FT
X_ph = angle(X);                        %Phase of the FT

figure(8);                              %Starting to plot
plot(f, X_ma, 'b'); axis tight;         %Starting to plot
title('Magnitude of the FT');           %Title
xlabel('f');                            %Title for the x-axis
print('images/exercice_4a','-dpng');    %Saving the plot

figure(9);                              %Starting to plot
plot(f, X_ph, 'c'); axis tight;         %Starting to plot
title('Phase of the FT');               %Title
xlabel('f');                            %Title for the x-axis
print('images/exercice_4b','-dpng');    %Saving the plot

%% MATCHED FILTER IN FREQUENCY DOMAIN %%

dt = 1e-5;                              %Sample spacing
df_i = 1;                               %Desired frequency resolution
t = -5e-3 : dt : 5e-3;                  %Time vector

u = 2*((t>=1e-3).*(t<=3e-3))-3*((t>=2e-3).*(t<=4e-3));	%Given signal
v = ((t>=-1e-3).*(t<=2e-3))+2*((t>=0).*(t<=1e-3));      %Given signal
s = u + 1i*v;                                           %Total signal

[X, f, df] = contFT(s,t(1),dt,df_i);    %Given funtion
X_ma = abs(X);                          %Magnitude of the FT
X_ph = angle(X);                        %Phase of the FT

figure(10);                             %Starting to plot
plot(f, X_ma, 'b'); axis tight;         %Starting to plot
title('Magnitude of the FT');           %Title
xlabel('f');                            %Title for the x-axis
print('images/exercice_5a','-dpng');    %Saving the plot

figure(11);                                 %Starting to plot
plot(f, X_ma, 'g'); xlim([-4e3,4e3]);       %Starting to plot
title('Magnitude of the FT (Zoomed-In)');   %Title
xlabel('f');                                %Title for the x-axis
print('images/exercice_5a_c','-dpng');      %Saving the plot

g = fliplr(conj(s));                    %Match filter
[h,y]= contconv(s,g,t(1),t(1),dt);      %Convolution
[H,k,~] = contFT(h,y(1),dt,df_i);       %Given funtion
H_ma = abs(H);                          %Magnitude of the FT
H_ph = angle(H);                        %Phase of the FT
X2 = abs(X).^2;                         %Magnitude square of X

figure(12);                                         %Starting to plot
subplot(2,1,1);                                     %First plot
plot(k, H_ma, 'g'); axis tight;                     %Plot
title('Magnitude of the FT of the convolution');    %Title
xlabel('f');                                        %Title for the x-axis
subplot(2,1,2);                                     %Second plot
plot(f, X2, 'c'); axis tight;                       %Plot
title('Magnitude square of the FT');                %Title
xlabel('f');                                        %Title for the x-axis
print('images/exercice_5b','-dpng');                %Saving the plot

figure(13);                                                 %Starting to plot
subplot(2,1,1);                                             %First plot
plot(k, H_ma, 'b'); xlim([-4e3,4e3]);                       %Plot
title('Magnitude of the FT of the convolution (Zoomed-In)');%Title
xlabel('f');                                                %Title for the x-axis
subplot(2,1,2);                                             %Second plot
plot(f, X2, 'c'); xlim([-4e3,4e3]);                         %Plot
title('Magnitude square of the FT (Zoomed-In)');            %Title
xlabel('f');                                                %Title for the x-axis
print('images/exercice_5b_c','-dpng');                      %Saving the plot

figure(14);                                         %Starting to plot
plot(k, H_ph, 'g'); axis tight;                     %Starting to plot
title('Phase of the FT of the convolution');        %Title
xlabel('f');                                        %Title for the x-axis
print('images/exercice_5c','-dpng');                %Saving the plot

figure(15);                                                 %Starting to plot
plot(k, H_ph, 'b'); ylim([-1,1]);                           %Starting to plot
title('Phase of the FT of the convolution (Zoomed-Out)');   %Title
xlabel('f');                                                %Title for the x-axis
print('images/exercice_5c_f','-dpng');                      %Saving the plot
