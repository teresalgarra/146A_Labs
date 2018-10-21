%% LAB 2.2: MODELING A LAMPPOST BASED BROADBAND NETWORK %%

%The goal of this lab is to illustrate how wireless multipath channels 
%can be modeled in complex baseband.

clear;
clc;
close all;

%% LABORATORY ASSIGNMENT %%

clear;
clc;
close all;

%%2.1%%

f = 60e9;                                       %Unlicensed spectrum at 5GHz
c = 3e8;                                        %Speed of light
n1 = 1;                                         %Refractive index of air
n2 = 8;                                         %Refractive index of ground
D = 200;                                        %Distance between lampposts
ht = 10;                                        %Height of the first lamppost
hr = 10;                                        %Height of the second lamppost
k = 2*pi*f/c;                                   %Wave Number

d1 = sqrt(D.^2 + (ht-hr).^2);                   %Trigonometric formula
d2 = sqrt(D.^2 + (ht+hr).^2);                   %Trigonometric formula

sin_i = ht/d1;                                  %Angle of incidence
phi_i = asin(sin_i);                            %Angle of incidence
sin_t = n1/n2*sin_i;                            %Angle of transmission
phi_t = asin(sin_t);                            %Angle of transmission
ci = cos(phi_i);                                %Cosine of the angle of incidence
ct = cos(phi_t);

A1 = 1;                                         %Amplitude of the direct path
A2 = (n1*ci-n2*ct)/(n1*ci+n2*ct);               %Amplitude of the indirect path

delay = d2/c - d1/c;                            %Delay spread
bc = 1/delay;                                   %Coherence Bandwidth

%%2.2%%

Dc = 100;                                       %Distance to the car
hc = 2;                                         %Height of the car

d1c = sqrt(Dc.^2 + (ht-hc).^2);                 %Trigonometric formula
d2c = sqrt(Dc.^2 + (ht+hc).^2);                 %Trigonometric formula

sin_ic = ht/d1c;                                %Angle of incidence
phi_ic = asin(sin_ic);                          %Angle of incidence
sin_tc = n1/n2*sin_ic;                          %Angle of transmission
phi_tc = asin(sin_tc);                          %Angle of transmission
cic = cos(phi_ic);                              %Cosine of the angle of incidence
ctc = cos(phi_tc);

ph1c = k*d1c;                                   %Phase of the direct path
ph2c = k*d2c;                                   %Phase of the reflected path
A1c = 1;                                        %Amplitude of the direct path
A2c = (n1*cic-n2*ctc)/(n1*cic+n2*ctc);          %Amplitude of the indirect path
hc = A1c*exp(1j*ph1c) + A2c*exp(1j*ph2c);       %Channel formula

delayc = d2c/c - d1c/c;                         %Delay spread
bcc = 1/delayc;                                 %Coherence Bandwidth

%%Fading and diversity for the backhaul link%%

%%2.3%%

dd = 0.0001;                                    %Step for the height vector
hr = 9.8:dd:10.2;                               %Height of the second lamppost

hnom1 = [];                                     %Nominal Channel
h1 = [];                                        %Complex Channel
for hh = 9.8:dd:10.2                            %Loop start
    d1 = sqrt(D.^2 + (ht-hh).^2);               %Trigonometric formula
    d2 = sqrt(D.^2 + (ht+hh).^2);               %Trigonometric formula
    ph1 = k*d1;                                 %Phase of the direct path
    ph2 = k*d2;                                 %Phase of the reflected path
    hnom1 = [hnom1, A1*exp(1j*ph1)];            %Nominal channel
    h1 = [h1, A1*exp(1j*ph1) + A2*exp(1j*ph2)];	%Complex Channel
end                                             %End of loop

npg1 = 20*log10(abs(h1)./abs(hnom1));           %Normalized Power Gain

figure;                                         %Starting to plot
plot(hr, npg1, 'r'); axis tight;                %Plot
title('Normalized Power Gain [dB]');            %Title
xlabel('h_{r}');                                %Time axis
print('images/2_3','-dpng');                    %Saving the plot

%%2.4%%

i = 0;                                          %Counter
for x = 1:length(npg1)                          %Going through npg1
    if npg1(x) <= -10                           %Condition
        i = i+1;                                %Adding counter
    end                                         %End of conditional loop
end                                             %End of for loop

probability = 100*i/length(npg1);               %Probability in %

%%2.5%%

ht1 = 10;                                       %Height of the first antenna 
ht2 = 10.1;                                     %Height of the second antenna 
hnom1 = [];                                     %Nominal Channel
hnom2 = [];                                     %Nominal Channel
h1 = [];
h2 = [];                                        %Complex Channel
for hh = 9.8:dd:10.2                            %Loop start
    d11 = sqrt(D.^2 + (ht1-hh).^2);             %Trigonometric formula
    d21 = sqrt(D.^2 + (ht1+hh).^2);             %Trigonometric formula
    ph11 = k*d11;                               %Phase of the direct path
    ph21 = k*d21;                               %Phase of the reflected path
    hnom1 = [hnom1, A1*exp(1j*ph11)];           %Nominal channel
    h1 = [h1,A1*exp(1j*ph11)+A2*exp(1j*ph21)];  %Complex Channel
    d12 = sqrt(D.^2 + (ht2-hh).^2);             %Trigonometric formula
    d22 = sqrt(D.^2 + (ht2+hh).^2);             %Trigonometric formula
    ph12 = k*d12;                               %Phase of the direct path
    ph22 = k*d22;                               %Phase of the reflected path
    hnom2 = [hnom2, A1*exp(1j*ph12)];           %Nominal channel
    h2 = [h2,A1*exp(1j*ph12)+A2*exp(1j*ph22)];  %Complex Channel
end                                             %End of loop

npg1 = 20*log10(abs(h1)./abs(hnom1));           %Normalized Power Gain
npg2 = 20*log10(abs(h2)./abs(hnom2));           %Normalized Power Gain

figure;                                         %Starting to plot
subplot(3,1,1);                                 %First plot
plot(hr, npg1, 'r'); axis tight;                %Plot
title('Normalized Power Gain [dB] for h1');     %Title
xlabel('h_{r}');                                %Time axis
subplot(3,1,2);                                 %Second plot
plot(hr, npg2, 'b'); axis tight;                %Plot
title('Normalized Power Gain [dB] for h2');     %Title
xlabel('h_{r}');                                %Time axis
subplot(3,1,3);                                 %Third plot
plot(hr, npg1, 'r'); axis tight;                %Plot
hold on;                                        %Holding the plot
plot(hr, npg2, 'b'); axis tight;                %Plot
hold off;                                       %Holding off
title('Normalized Power Gain [dB] for h');      %Title
xlabel('h_{r}');                                %Time axis
print('images/2_5','-dpng');                    %Saving the plot

%%2.6%%

h3 = max(h1, h2);                               %Maximum of h
hnom3 = max(hnom1, hnom2);                      %Maximum of hnom
npg3 = 20*log10(abs(h3)./abs(hnom2));           %Normalized Power Gain

figure;                                         %Starting to plot
plot(hr, npg3, 'r'); axis tight;                %Plot
title('Normalized Power Gain [dB]');            %Title
xlabel('h_{r}');                                %Time axis

%%2.7%%

i = 0;                                          %Counter
for x = 1:length(npg3)                          %Going through npg3
    if npg3(x) <= -10                           %Condition
        i = i+1;                                %Adding counter
    end                                         %End of conditional loop
end                                             %End of for loop

probability2 = 100*i/length(npg3);              %Probability in %
