dt = 0.001;
bc = 0.47619;
t = -2*bc:dt:2*bc;

h = 2.*exp(-1j.*2.*pi.*t.*0.1)+1j.*exp(-1j.*2.*pi.*t.*0.64)-0.8.*exp(-1j.*2.*pi.*t.*2.2);

h_ma = abs(h);
h_ph = angle(h);

grid on;

figure;
subplot(2,1,1);
plot(t, h_ma, 'c'); xlim([-0.75, -0.25]);
title('Magnitude of H(f)');
subplot(2,1,2);
plot(t, h_ph, 'g'); xlim([-0.75, -0.25]);
title('Phase of H(f)');
print('13_close','-dpng');

figure;
subplot(2,1,1);
plot(t, h_ma, 'c');
title('Magnitude of H(f)');
subplot(2,1,2);
plot(t, h_ph, 'g');
title('Phase of H(f)');
print('13','-dpng');

y = 10*log10(h_ma);

figure;
plot(t,y);
title('Magnitude of H(f) in dB');
print('13_c','-dpng');

sum = 0;

for n = -2*bc/dt : dt : 2*bc/dt
    sum = sum + abs(h).^2.*dt;
end

figure;
plot(t,sum);
title('Average Channel Power');
print('13_c2','-dpng');

sum2 = 0;

for n = -2/dt : dt : 2/dt
    sum2 = sum2 + abs(h).^2.*dt;
end

figure;
plot(t,sum2);
title('Average Channel Power Normalized');
print('13_d','-dpng');