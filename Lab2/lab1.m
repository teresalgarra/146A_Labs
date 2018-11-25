%% LAB 2.1: MODELING CARRIER PHASE UNCERTAINTY %%

%The goal of this lab is to explore modeling and receiver operations in 
%complex baseband. In particular, we model and undo the effect of carrier 
%phase mismatch between the receiver LO and the incoming carrier.

clear;
clc;
close all;

%% LABORATORY ASSIGNMENT %%

%%Random coeficients%%
set = [-1,1];                           %These are the coefficients
dt = 0.01;                              %Step for discrete time
osfactor = 1/dt;                        %Over Sampling factor
N = 100;                                %Samples
tt = dt:dt:N;                           %Time vector

%%1.1: uc and us%%

u_c = [];                               %Pre-allocating
for n = 1:1:N                           %Starting the loop
    pos_c = randi(length(set));         %Random position of the set for bc
    b_c = set(pos_c);                   %Coefficient bc
    bcc = repmat(b_c, osfactor, 1);     %Matrix of bc
    bbc = bcc(:)';                      %Changing size
    u_c = [u_c, bbc];                   %Putting it in uc
end

u_s = [];                               %Pre-allocating
for n = 1:1:N                           %Starting the loop
    pos_s = randi(length(set));         %Random position of the set for bs
    b_s = set(pos_s);                   %Coefficient bs
    bss = repmat(b_s, osfactor, 1);     %Matrix of bs
    bbs = bss(:)';                      %Changing size
    u_s = [u_s, bbs];                   %Putting it in us
end

figure;                                 %Starting to plot
subplot(2,1,1);                         %First plot
plot(tt, u_c, 'b'); xlim([1,10]);       %Plot
title('Uc');                            %Title
xlabel('t');                            %Time axis
subplot(2,1,2);                         %Second plot
plot(tt, u_s, 'r'); xlim([1,10]);       %Plot
title('Us');                            %Title
xlabel('t');                            %Time axis
print('images/1_1','-dpng');            %Saving the plot

%%1.2: Upconversion%%


u_p1 = u_c.*cos(40.*pi.*tt);            %Given formula

figure;                                 %Starting to plot
plot(tt, u_p1, 'b'); xlim([1,4]);       %Plot
title('Up without Q');                  %Title
xlabel('t');                            %Time axis
print('images/1_2','-dpng');            %Saving the plot

%%1.3: Passband signal%%

u_p = u_p1 - u_s.*sin(40*pi*tt);        %Given formula

figure;                                 %Starting to plot
plot(tt, u_p, 'r'); xlim([1,4]);        %Plot
title('Up with Q');                     %Title
xlabel('t');                            %Time axis
print('images/1_3','-dpng');            %Saving the plot

%%1.4: Downconversion%%

lpf = ((tt>=0).*(tt<=0.25));                %Low Pass Filter
v_c1 = u_p.*2.*cos(40*pi*tt);               %Given formula
[v_c,t_c] = contconv(v_c1,lpf,0,0,0.01);    %Convolution
v_s1 = u_p.*2.*sin(40*pi*tt);               %Given formula
[v_s,t_s] = contconv(v_s1,lpf,0,0,0.01);    %Convolution

figure;                                 %Starting to plot
subplot(2,1,1);                         %First plot
plot(t_c, v_c, 'b'); xlim([1,N]);       %Plot
title('Vc');                            %Title
xlabel('t');                            %Time axis
subplot(2,1,2);                         %Second plot
plot(t_s, v_s, 'r'); xlim([1,N]);       %Plot
title('Vs');                            %Title
xlabel('t');                            %Time axis
print('images/1_4','-dpng');            %Saving the plot

%%1.5: Downconversion with phase%%

w_c1 = u_p.*2.*cos(40*pi*tt + pi/4);        %Given formula
[w_c,tt_c] = contconv(w_c1,lpf,0,0,0.01);   %Convolution
w_s1 = u_p.*2.*sin(40*pi*tt + pi/4);        %Given formula
[w_s,tt_s] = contconv(w_s1,lpf,0,0,0.01);   %Convolution

figure;                                 %Starting to plot
subplot(2,1,1);                         %First plot
plot(tt_c, w_c, 'b'); xlim([1,N]);      %Plot
title('Vc');                            %Title
xlabel('t');                            %Time axis
subplot(2,1,2);                         %Second plot
plot(tt_s, w_s, 'r'); xlim([1,N]);      %Plot
title('Vs');                            %Title
xlabel('t');                            %Time axis
print('images/1_5','-dpng');            %Saving the plot

%%1.6: Recovering uc and us%%

u_c2 = real(w_c * exp(1j*pi/4));        %Exponential needed to recover uc
u_s2 = real(w_s * exp(1j*pi/4));        %Exponential needed to recover us

figure;                                 %Starting to plot
subplot(2,1,1);                         %First plot
plot(tt_c, u_c2, 'b'); xlim([1,N]);     %Plot
title('Recovered Uc');                  %Title
xlabel('t');                            %Time axis
subplot(2,1,2);                         %Second plot
plot(tt_s, u_s2, 'r'); xlim([1,N]);     %Plot
title('Recovered Us');                  %Title
xlabel('t');                            %Time axis
print('images/1_6','-dpng');            %Saving the plot

figure;                                 %Starting to plot
subplot(2,1,1);                         %First plot
plot(tt_c, u_c2, 'b'); xlim([1,N]);     %Plot
title('Recovered Uc');                  %Title
xlabel('t');                            %Time axis
subplot(2,1,2);                         %Second plot
plot(tt, u_c, 'r'); xlim([1,N]);        %Plot
title('Uc');                            %Title
xlabel('t');                            %Time axis
print('images/1_6a','-dpng');           %Saving the plot

figure;                                 %Starting to plot
subplot(2,1,1);                         %First plot
plot(tt_c, u_s2, 'b'); xlim([1,N]);     %Plot
title('Recovered Us');                  %Title
xlabel('t');                            %Time axis
subplot(2,1,2);                         %Second plot
plot(tt, u_s, 'r'); xlim([1,N]);        %Plot
title('Us');                            %Title
xlabel('t');                            %Time axis
print('images/1_6b','-dpng');           %Saving the plot

