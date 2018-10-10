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
N=100;                                %Number of samples
tt= 1:1:N;                            %Time vector for plotting
p = linspace(0,N);                    %Unitary pulse
p(1:2) = 1;                           %It's 1 until x=1
p(3:N) = 0;                           %Else, it's 0.

%%Signals%%

symbols = 10;

%u_c%
u_c = zeros(N, 1);
for i = 1:1:symbols
    for n = 1:1:N
        for t = 1:1:N
            if t>n
                pos_c = randi(length(set));           %This selects a random position of the set for b_c
                b_c = set(pos_c);                     %We get the coefficient b_c
                u_c(t) = b_c.*p(t-n);                 %We create the function u_c
            end
        end
    end
end

%u_s%
u_s = zeros(N,1);
for i = 1:1:symbols
    for n = 1:1:N
        for t = 1:1:N
            if t>n
                pos_s = randi(length(set));           %This selects a random position of the set for b_c
                b_s = set(pos_s);                     %We get the coefficient b_s
                u_s(t) = b_s.*p(t-n);                 %We create the function u_s
            end
        end
    end
end

%%Plots%%
figure;
subplot(2,1,1);
plot(tt, u_c, 'b'); axis tight;
title('Uc');
xlabel('t');
subplot(2,1,2);
plot(tt, u_s, 'r'); axis tight;
title('Us');
xlabel('t');
print('~/Documents/ECE146A/Lab1/uc_us','-dpng');

%%Upconversion%%
symbols = 4;

%u_p without Q%
u_p = zeros(N, 1);
u_c = zeros(N, 1);
for i = 1:1:symbols
    for n = 1:1:N
        for t = 1:1:N
            if t>n
                pos_c = randi(length(set));           %This selects a random position of the set for b_c
                b_c = set(pos_c);                     %We get the coefficient b_c
                u_c(t) = b_c.*p(t-n);                 %We create the function u_c
                u_p(t) = u_c(t).*cos(40*pi*t);       %We apply the formula
            end
        end
    end
end

%%Passband signal%%

%u_p with Q%
u_p_t = zeros(N, 1);
u_c = zeros(N, 1);
u_s = zeros(N, 1);
for i = 1:1:symbols
    for n = 1:1:N
        for t = 1:1:N
            if t>n
                pos_c = randi(length(set));           %This selects a random position of the set for b_c
                b_c = set(pos_c);                     %We get the coefficient b_c
                u_c(t) = b_c.*p(t-n);                 %We create the function u_c
                pos_s = randi(length(set));           %This selects a random position of the set for b_c
                b_s = set(pos_s);                     %We get the coefficient b_s
                u_s(t) =  b_s.*p(t-n);                %We create the function u_s
                u_p_t(t) = u_c(t).*cos(40*pi*t) - u_s(t).*sin(40*pi*t);   %We apply the formula
            end
        end
    end
end

%%Plots%%
figure;
subplot(2,1,1);
plot(tt, u_p, 'b'); axis tight;
title('Up without Q');
xlabel('t');
subplot(2,1,2);
plot(tt, u_p_t, 'r'); axis tight;
title('Up with Q');
xlabel('t');
print('~/Documents/ECE146A/Lab1/up','-dpng');

