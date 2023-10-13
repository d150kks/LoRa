clear all
close all
clc

%% ================================= Variables
fc = 1000e6;
fs = 30e6;

M = 4;   
Pream_len = 4096;
Chip_len = 295;
Data_len = 12000;

% Формируем преамбулу
x = randi([0, M-1],1,Pream_len); % Случайное сообщение
% Data
data = randi([0, 1],1,Data_len); % Случайное сообщение


%% ================================= Modulation
% Используем 16-точечную КАМ
y = qammod(x, M);

pream = [y, y];

% формируем сигналы для 0 и 1
chip_zero = (2*randi([0,1],1,Chip_len)-1) + 1i*(2*randi([0,1],1,Chip_len)-1);
chip_one = (2*randi([0,1],1,Chip_len)-1) + 1i*(2*randi([0,1],1,Chip_len)-1);


sig_mod = zeros(1,Data_len*Chip_len);
for i = 1:Data_len
    if(data(i) == 1)
        sig_mod(i*Chip_len-Chip_len+1:i*Chip_len) = chip_one;
    end
    
    if(data(i) == 0)
        sig_mod(i*Chip_len-Chip_len+1:i*Chip_len) = chip_zero;
    end
end


frame_tx = [pream, sig_mod];

%% ================================= Bits per Second
Tsig = length(sig_mod)/fs;
bps = Data_len/Tsig;

fprintf('T = %.2f\n', Tsig)
fprintf('Bits per Sec = %.2f\n', bps)

figure(1)
plot(real(frame_tx))
return
