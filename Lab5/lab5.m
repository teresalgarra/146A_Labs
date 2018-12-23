%% LAB 5: AN INITIAL EXPOSURE TO MM-WAVE RADAR SENSING %%

%The goal of this lab is to provide an initial exposure to the fundamentals of
%radar sensing, with parameters consistent with emerging mm wave sensing
%applications.

clear;
clc;
close all;

%% SNR = 1 %%

clear;
clc;
close all;

S = 100e12;                                 %Slew rate
Td = 10e-6;                                 %Chirp duration
R1 = 3;                                     %Distance to target 1
R2 = 10;                                    %Distance to target 2
R3 = 30;                                    %Distance to target 3
c = 3*1e8;                                  %Speed of light
tau1 = 2*R1/c;                              %Delay of first target
tau2 = 2*R2/c;                              %Delay of second target
tau3 = 2*R3/c;                              %Delay of third target
f1 = S*tau1;                                %Frequency of first target
f2 = S*tau2;                                %Frequency of second target
f3 = S*tau3;                                %Frequency of third target
fs = 2.5*f3;                                %Frequency of sampling
Ts = 1/fs;                                  %Period of sampling
t = 0:Ts:Td;                                %Time vector
N = length(t);                              %Number of samples
SNR1 = 1;                                   %SNR to target 1
SNR2 = SNR1*(R1/R2).^4;                     %SNR to target 2
SNR3 = SNR1*(R1/R3).^4;                     %SNR to target 3
A1 = sqrt(SNR1/N);                          %Amplitude of target 1
A2 = sqrt(SNR2/N);                          %Amplitude of target 2
A3 = sqrt(SNR3/N);                          %Amplitude of target 3
nfft = nextpow2(N)*64;                      %Size of the Fourier Transform
fi = 0:1/(nfft*Ts):(nfft-1)/(nfft*Ts);      %Frequency axis
iterations = 100;                           %Number of iterations

e1 = [];                                    %Preallocating error vector 1
e2 = [];                                    %Preallocating error vector 2
e3 = [];                                    %Preallocating error vector 3

for q = 1:iterations
    
    rng('shuffle');                         %Changing the seed for rand

    y1 = A1.*exp(1j*2*pi.*f1.*t);           %First target
    y2 = A2.*exp(1j*2*pi.*f2.*t);           %Second target
    y3 = A3.*exp(1j*2*pi.*f3.*t);           %Third target
    n = 0.5.*(randn(1,N)+1j.*randn(1,N));   %Noise

    x = y1 + y2 + y3;                       %No-noise signal
    X = fft(x, nfft);                       %Fourier Transform
    X = abs(X);                             %Absolute Value
    y = y1 + y2 + y3 + n;                   %Total signal
    Y = fft(y, nfft);                       %Fourier transform
    Y = abs(Y);                             %Absolute Value

    [peaks, peaks_i] = find_peaks(Y);       %Find Peaks function

    for h = 1:length(peaks)-1               %First sorting loop
        for i = 1:length(peaks)-1           %Second sorting loop
            if peaks(i)< peaks(i+1)         %Condition for sorting
                comodin = peaks(i);         %Sorting
                peaks(i) = peaks(i+1);      %Sorting
                peaks(i+1) = comodin;       %Sorting
                comodin_i = peaks_i(i);     %Sorting
                peaks_i(i) = peaks_i(i+1);  %Sorting
                peaks_i(i+1) = comodin_i;   %Sorting
            end                             %End of contidion   
        end                                 %End of second loop
    end                                     %End of first loop

    peaks_5 = peaks(:, 1:5);                %5 highest peaks
    peaks_i_5 = peaks_i(:, 1:5);            %Location of 5 highest peaks

    index_5 = zeros(1,5);                   %Index vector
    for i = 1:5                             %Loop
        index_5(i) = fi(peaks_i_5(1,i));    %Frequency of every peak location
    end                                     %End of loop
    
    range = index_5./S.*c./2;               %Range estimate

    e1 = [e1 min(abs(range - R1))];         %Error for target 1
    e2 = [e2 min(abs(range - R2))];         %Error for target 2
    e3 = [e3 min(abs(range - R3))];         %Error for target 3

end

figure;                                     %Plotting
subplot(2,1,1);                             %First plot
plot(fi, X, 'b');                           %No-noise signal
xlabel('f');                                %X-axis
title('No-noise Signal');                   %Title
subplot(2,1,2);                             %Second plot
plot(fi, Y, 'r');                           %Total signal
xlabel('f');                                %X-axis
title('Signal with noise');                 %Title
print('images/signals_1','-dpng');          %Saving the plot

L = 100;                                    %Beans for the histogram

figure;                                     %Plotting
hist(e1, L);                                %Histogram for error1
xlabel('Range');                            %X-axis
title('Error for Target 1');                %Title
print('images/hist_e1_1','-dpng');          %Saving the plot

figure;                                     %Plotting
hist(e2, L);                                %Histogram for error2
xlabel('Range');                            %X-axis
title('Error for Target 2');                %Title
print('images/hist_e2_1','-dpng');          %Saving the plot

figure;                                     %Plotting
hist(e3, L);                                %Histogram for error3
xlabel('Range');                            %X-axis
title('Error for Target 3');                %Title
print('images/hist_e3_1','-dpng');          %Saving the plot

%% SNR = 10 %%

clear;
clc;
close all;

fc = 80e9;                                  %Reference frequency
S = 100e12;                                 %Slew rate
Td = 10e-6;                                 %Chirp duration
R1 = 3;                                     %Distance to target 1
R2 = 10;                                    %Distance to target 2
R3 = 30;                                    %Distance to target 3
c = 3*1e8;                                  %Speed of light
tau1 = 2*R1/c;                              %Delay of first target
tau2 = 2*R2/c;                              %Delay of second target
tau3 = 2*R3/c;                              %Delay of third target
f1 = S*tau1;                                %Frequency of first target
f2 = S*tau2;                                %Frequency of second target
f3 = S*tau3;                                %Frequency of third target
fs = 2.5*f3;                                %Frequency of sampling
Ts = 1/fs;                                  %Period of sampling
t = 0:Ts:Td;                                %Time vector
N = length(t);                              %Number of samples
SNR1 = 10;                                  %SNR to target 1
SNR2 = SNR1*(R1/R2).^4;                     %SNR to target 2
SNR3 = SNR1*(R1/R3).^4;                     %SNR to target 3
A1 = sqrt(SNR1/N);                          %Amplitude of target 1
A2 = sqrt(SNR2/N);                          %Amplitude of target 2
A3 = sqrt(SNR3/N);                          %Amplitude of target 3
nfft = nextpow2(N)*64;                      %Size of the Fourier Transform
fi = 0:1/(nfft*Ts):(nfft-1)/(nfft*Ts);      %Frequency axis
iterations = 100;                           %Number of iterations

e1 = [];                                    %Preallocating error vector 1
e2 = [];                                    %Preallocating error vector 2
e3 = [];                                    %Preallocating error vector 3

for q = 1:iterations
    
    rng('shuffle');                         %Changing the seed for rand

    y1 = A1.*exp(1j*2*pi.*f1.*t);           %First target
    y2 = A2.*exp(1j*2*pi.*f2.*t);           %Second target
    y3 = A3.*exp(1j*2*pi.*f3.*t);           %Third target
    n = 0.5.*(randn(1,N)+1j.*randn(1,N));   %Noise

    x = y1 + y2 + y3;                       %No-noise signal
    X = fft(x, nfft);                       %Fourier Transform
    X = abs(X);                             %Absolute Value
    y = y1 + y2 + y3 + n;                   %Total signal
    Y = fft(y, nfft);                       %Fourier transform
    Y = abs(Y);                             %Absolute Value

    [peaks, peaks_i] = find_peaks(Y);       %Find Peaks function

    for h = 1:length(peaks)-1               %First sorting loop
        for i = 1:length(peaks)-1           %Second sorting loop
            if peaks(i)< peaks(i+1)         %Condition for sorting
                comodin = peaks(i);         %Sorting
                peaks(i) = peaks(i+1);      %Sorting
                peaks(i+1) = comodin;       %Sorting
                comodin_i = peaks_i(i);     %Sorting
                peaks_i(i) = peaks_i(i+1);  %Sorting
                peaks_i(i+1) = comodin_i;   %Sorting
            end                             %End of contidion  
        end                                 %End of second loop
    end                                     %End of first loop

    peaks_5 = peaks(:, 1:5);                %5 highest peaks
    peaks_i_5 = peaks_i(:, 1:5);            %Location of 5 highest peaks

    index_5 = zeros(1,5);                   %Index vector
    for i = 1:5                             %Loop
        index_5(i) = fi(peaks_i_5(1,i));    %Frequency of every peak location
    end                                     %End of loop
    
    range = index_5./S.*c./2;               %Range estimate

    e1 = [e1 min(abs(range - R1))];         %Error for target 1
    e2 = [e2 min(abs(range - R2))];         %Error for target 2
    e3 = [e3 min(abs(range - R3))];         %Error for target 3

end

figure;                                     %Plotting
subplot(2,1,1);                             %First plot
plot(fi, X, 'b');                           %No-noise signal
xlabel('f');                                %X-axis
title('No-noise Signal');                   %Title
subplot(2,1,2);                             %Second plot
plot(fi, Y, 'r');                           %Total signal
xlabel('f');                                %X-axis
title('Signal with noise');                 %Title
print('images/signals_10','-dpng');         %Saving the plot

L = 100;                                    %Beans for the histogram

figure;                                     %Plotting
hist(e1, L);                                %Histogram for error1
xlabel('Range');                            %X-axis
title('Error for Target 1');                %Title
print('images/hist_e1_10','-dpng');         %Saving the plot

figure;                                     %Plotting
hist(e2, L);                                %Histogram for error2
xlabel('Range');                            %X-axis
title('Error for Target 2');                %Title
print('images/hist_e2_10','-dpng');         %Saving the plot

figure;                                     %Plotting
hist(e3, L);                                %Histogram for error3
xlabel('Range');                            %X-axis
title('Error for Target 3');                %Title
print('images/hist_e3_10','-dpng');         %Saving the plot

%% SNR = 50 %%

clear;
clc;
close all;

fc = 80e9;                                  %Reference frequency
S = 100e12;                                 %Slew rate
Td = 10e-6;                                 %Chirp duration
R1 = 3;                                     %Distance to target 1
R2 = 10;                                    %Distance to target 2
R3 = 30;                                    %Distance to target 3
c = 3*1e8;                                  %Speed of light
tau1 = 2*R1/c;                              %Delay of first target
tau2 = 2*R2/c;                              %Delay of second target
tau3 = 2*R3/c;                              %Delay of third target
f1 = S*tau1;                                %Frequency of first target
f2 = S*tau2;                                %Frequency of second target
f3 = S*tau3;                                %Frequency of third target
fs = 2.5*f3;                                %Frequency of sampling
Ts = 1/fs;                                  %Period of sampling
t = 0:Ts:Td;                                %Time vector
N = length(t);                              %Number of samples
SNR1 = 50;                                  %SNR to target 1
SNR2 = SNR1*(R1/R2).^4;                     %SNR to target 2
SNR3 = SNR1*(R1/R3).^4;                     %SNR to target 3
A1 = sqrt(SNR1/N);                          %Amplitude of target 1
A2 = sqrt(SNR2/N);                          %Amplitude of target 2
A3 = sqrt(SNR3/N);                          %Amplitude of target 3
nfft = nextpow2(N)*64;                      %Size of the Fourier Transform
fi = 0:1/(nfft*Ts):(nfft-1)/(nfft*Ts);      %Frequency axis
iterations = 100;                           %Number of iterations

e1 = [];                                    %Preallocating error vector 1
e2 = [];                                    %Preallocating error vector 2
e3 = [];                                    %Preallocating error vector 3

for q = 1:iterations
    
    rng('shuffle');                         %Changing the seed for rand

    y1 = A1.*exp(1j*2*pi.*f1.*t);           %First target
    y2 = A2.*exp(1j*2*pi.*f2.*t);           %Second target
    y3 = A3.*exp(1j*2*pi.*f3.*t);           %Third target
    n = 0.5.*(randn(1,N)+1j.*randn(1,N));   %Noise

    x = y1 + y2 + y3;                       %No-noise signal
    X = fft(x, nfft);                       %Fourier Transform
    X = abs(X);                             %Absolute Value
    y = y1 + y2 + y3 + n;                   %Total signal
    Y = fft(y, nfft);                       %Fourier transform
    Y = abs(Y);                             %Absolute Value

    [peaks, peaks_i] = find_peaks(Y);       %Find Peaks function

    for h = 1:length(peaks)-1               %First sorting loop
        for i = 1:length(peaks)-1           %Second sorting loop
            if peaks(i)< peaks(i+1)         %Condition for sorting
                comodin = peaks(i);         %Sorting
                peaks(i) = peaks(i+1);      %Sorting
                peaks(i+1) = comodin;       %Sorting
                comodin_i = peaks_i(i);     %Sorting
                peaks_i(i) = peaks_i(i+1);  %Sorting
                peaks_i(i+1) = comodin_i;   %Sorting
            end                             %End of contidion  
        end                                 %End of second loop
    end                                     %End of first loop

    peaks_5 = peaks(:, 1:5);                %5 highest peaks
    peaks_i_5 = peaks_i(:, 1:5);            %Location of 5 highest peaks

    index_5 = zeros(1,5);                   %Index vector
    for i = 1:5                             %Loop
        index_5(i) = fi(peaks_i_5(1,i));    %Frequency of every peak location
    end                                     %End of loop
    
    range = index_5./S.*c./2;               %Range estimate

    e1 = [e1 min(abs(range - R1))];         %Error for target 1
    e2 = [e2 min(abs(range - R2))];         %Error for target 2
    e3 = [e3 min(abs(range - R3))];         %Error for target 3

end

figure;                                     %Plotting
subplot(2,1,1);                             %First plot
plot(fi, X, 'b');                           %No-noise signal
xlabel('f');                                %X-axis
title('No-noise Signal');                   %Title
subplot(2,1,2);                             %Second plot
plot(fi, Y, 'r');                           %Total signal
xlabel('f');                                %X-axis
title('Signal with noise');                 %Title
print('images/signals_50','-dpng');         %Saving the plot

L = 100;                                    %Beans for the histogram

figure;                                     %Plotting
hist(e1, L);                                %Histogram for error1
xlabel('Range');                            %X-axis
title('Error for Target 1');                %Title
print('images/hist_e1_50','-dpng');         %Saving the plot

figure;                                     %Plotting
hist(e2, L);                                %Histogram for error2
xlabel('Range');                            %X-axis
title('Error for Target 2');                %Title
print('images/hist_e2_50','-dpng');         %Saving the plot

figure;                                     %Plotting
hist(e3, L);                                %Histogram for error3
xlabel('Range');                            %X-axis
title('Error for Target 3');                %Title
print('images/hist_e3_50','-dpng');         %Saving the plot