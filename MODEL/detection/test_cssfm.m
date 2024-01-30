clc
clear all
close all

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 2e6;
snr = [-20:1:0];
nIter = 20;

LORA = myLoRaClass_true(SF,BW); 
LORA.fir_win = 1;
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;
Ts = LORA.Ts;

num_pre = 4;
intrlv_state = 13;

% num_sym = 1000; 
% numinfobits = num_sym*bits2sym; 
% data = randi([0 1],1, numinfobits); 

%% ================================= LDPC coding
numinfobits = 7;

% Message/Iformation bits
data = randi([0 1],1, numinfobits); 
num_sym = length(data)/SF;

%% ================================= Mодуляция
% RS
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data, num_sym, 1);

tx_preamble = repmat(chirp,1,num_pre);
tx_downch = repmat(downch,1,num_pre);

tx_chirp = [tx_downch, tx_preamble, mod_chirp];
tx_length = length(tx_chirp);



% SF = 7; 
% Base = 2^SF; 
% BW = 125e3; 
% Ts = (2^SF)/BW;
% ts = 1/(BW);              % время дискретизации
% m = BW/Ts;              % коэффициент измнения частоты
% t = 0:ts:Ts-ts;         % полоса времени ОДНОГО чирпа
% 
% % ~~~~~~~~ Ref Chirps ~~~~~~~~
% chirp  = exp( 1i * (2*pi*(m*t.^2)/2) );
% downch = exp(-1i * (2*pi*(m*t.^2)/2) );

%% ================================= Канал
sig_rx = awgn(mod_chirp, 100, 'measured');


I_Rx = real(sig_rx);     % Синфазная составляющая комплексной огибающей
Q_Rx = imag(sig_rx);     % Квадратурная составляющая комплексной огибающей
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% ~~~~~~~~ Demodulation the signal ~~~~~~~~
dI_Rx = [0, diff(I_Rx)*1]; % Вычисление дифференциалов
dQ_Rx = [0, diff(Q_Rx)*1]; % Вычисление дифференциалов

W_t_demod = (dQ_Rx.*I_Rx - Q_Rx.*dI_Rx)./(I_Rx.^2 + Q_Rx.^2);
Sm_demod = W_t_demod/1;


figure(1); hold on
plot(abs(fft(sig_rx(1:128).*downch)))
% plot(Sm_demod)
plot(abs(fft(Sm_demod)))











