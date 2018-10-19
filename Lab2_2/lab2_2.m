%% LAB 2.2: MODELING A LAMPPOST BASED BROADBAND NETWORK %%

%The goal of this lab is to illustrate how wireless multipath channels 
%can be modeled in complex baseband.

clear;
clc;
close all;

%% LABORATORY ASSIGNMENT %%

f = 60e9;                       %Unlicensed spectrum at 5GHz
A1 = 1;                         %Direct path
phi_1 = 0;                        %Direct path

n1 = 1;                         %Refractive index of air
n2 = 8;                         %Refractive index of ground

D = 200;
ht = 10;
hr = 10;

d1 = sqrt(D.^2 + (ht-hr).^2);
d2 = sqrt(D.^2 + (ht+hr).^2);

sin_i = ht/d1;
phi_i = asin(sin_i);
sin_t = n1/n2*sin_i;
phi_t = asin(sin_t);

sin_3 = (ht+hr)/d2;
phi_3 = asin(sin_3);
phi_2 = pi/2 - phi_3;


A2 = (n1*cos(phi_i)-n2*cos(phi_t))/(n1*cos(phi_i)+n2*cos(phi_t));

h = A1*exp(1j*phi_1) + A2*exp(1j*phi_2);

%%2.1%%

delay = phi_2 - phi_1;
bc = 1/delay;

%%2.2%%

%The delay spread 

%% Fading and diversity for the backhaul link %%