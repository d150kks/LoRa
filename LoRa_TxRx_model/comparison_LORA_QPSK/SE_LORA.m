clc
clear all
close all


tic
%% ================================= Переменные
% коэффициенты
SF = 10;        % коэффициент расширения спектра (от 7 до 12)
Base = 2^SF;   % База сигнала
num_sym = 1000;% Число передаваемых символов
num_pre = 4;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации

% массивы данных
numbits = (SF-4)*num_sym;          % число бит
data = randi([0, 1], 1, numbits); % формирование массива бит

% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 20e6+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 1000e6;            % центральная частота
band = BW;
Ts = (2^SF)/BW;        % длительность сигнала


%% ================================= модуляция
os = 1;
fs = BW*os;
ts=1/fs;

[mod_chirp, check_chirp, downch, check_data] = lorax_modified( data, fs, Base, Ts, SF, BW);
N = length(downch);

% rms(real(mod_chirp))
% mod_chirp = qammod(data.',4,'InputType','bit');
% peak = max(real(mod_chirp).^2);
% avr = mean(real(mod_chirp).^2);
% PAPR = 10*log10(peak/avr)
% return
%% Forming Signal
preamble = repmat(check_chirp,1,num_pre);
%%tx_chirp = [preamble, mod_chirp];
tx_chirp = [downch, preamble, mod_chirp];
tx_length = length(tx_chirp);


%% ================================= Bits per Second
CR = 3/4;
Tsig = Ts*num_sym;
bps = CR*numbits/Tsig;

fprintf('T = %.2f\n', Tsig)
fprintf('Bits per Sec = %.2f\n', bps)
% SF * (BW/2^SF) 

