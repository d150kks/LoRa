clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
SF = 10;        % коэффициент расширения спектра (от 7 до 12)
Base = 2^SF;   % База сигнала
num_sym = 10;% Число передаваемых символов
num_pre = 8;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
snr = 10;      % ОСШ

% массивы данных
numbits = 486;
data = randi([0, 1], 1, numbits);

% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 125e3+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 2.4e9;            % центральная частота
Ts = (Base)/BW;        % длительность сигнала


%% ================================= Mодуляция
os = 1;
fs = BW*os;
ts=1/fs;

[mod_chirp, check_chirp, downch, check_data] = lorax_modified( data, fs, Base, Ts, SF, BW);
N = length(downch);


%% ================================= Polar coding
E = 512;
numbits = 384;
data = randi([0, 1], 1, numbits);
codeword = nrPolarEncode(data.',E,10,false);
% codeword = nrPolarEncode(data.',E);


codeword_bad=codeword;


% nVar = 0.0; 
% chan = comm.AWGNChannel('NoiseMethod','Variance','Variance',nVar);
% mod = nrSymbolModulate(codeword_bad,'QPSK');
% rSig = chan(mod);
% 
% rxLLR = nrSymbolDemodulate(rSig,'QPSK',nVar);
% return

% 
codeword_bad(10) = ~codeword_bad(10);
codeword_bad(131) = ~codeword_bad(131);
codeword_bad(256) = ~codeword_bad(256);
codeword_bad(366) = ~codeword_bad(366);
codeword_bad(406) = ~codeword_bad(406);
codeword_bad = 2*~codeword_bad-1;

% codeword_bad(codeword_bad>0)=20;
% codeword_bad(codeword_bad<= 0)=-20;
%% ================================= LDPC decoding
L = 64;
rxbits = nrPolarDecode(codeword_bad,numbits,E,L,10,false,11).';
% L = 8;
% rxbits = nrPolarDecode(codeword_bad,numbits,E,L).';
% rxbits = nrPolarDecode(rxLLR,numbits,E,L).';

errStats = sum(data~=rxbits);
fprintf('Number of errors = %d\n', errStats)


figure(1)
plot(normalize(2*~codeword-1), 'LineWidth',2)
hold on
% plot(normalize(codeword_bad))