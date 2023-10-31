clc
clear all
close all




%% ================================= Переменные
ts = pi/200;
t = 0:ts:2*pi-ts;
s = cos(t*20);
len = length(s);
s_rect = s.*fftshift(gausswin(len).');

figure(1)
subplot(221)
plot(s)

subplot(222)
stem(abs(fft(s)))

subplot(223)
plot(s_rect)

subplot(224)
plot(abs(fft(s_rect)))

