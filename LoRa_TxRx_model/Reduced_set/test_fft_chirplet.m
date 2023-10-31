clc
clear all
close all


%% ================================= Переменные
N = 128;
y = zeros(1,N);
% x;
% for k=1:N
%     for j=1:N
%         y(k) = y(k)+

j=1:N;
k = 128;
m=1/(N^2);
d=exp(-1i*2*pi*k*j/N);
a=exp( 1i*2*pi*k*m*j.^2/2 );

figure(1); hold on
plot(real(d))
plot(real(a))