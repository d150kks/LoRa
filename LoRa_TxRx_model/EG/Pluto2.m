clear all
close all
clc

fc = 1000e6;
fs = 30e6;

M = 4;   
Pream_len = 1000;
Chip_len = 400;
Data_len = 400;

% Формируем преамбулу
x = randi([0, M-1],1,Pream_len); % Случайное сообщение

% Используем 16-точечную КАМ
y = qammod(x, M);
y(1:100)=0;
y(901:end)=0;
y = fftshift(y);
y(1)=0;

Buf_len_koef = 3;
figure
plot(abs(y));

sig1 = ifft(y);
pream = [sig1, sig1];
% формируем сигналы для 0 и 1


% zero chip

chip_zero = randn(1,Chip_len) + 1i*randn(1,Chip_len);
chip_zero = chip_zero/max(chip_zero);

% one chip

chip_one = randn(1,Chip_len) + 1i*randn(1,Chip_len);
chip_one = chip_one/max(chip_one);

data = randi([0, 1],1,Data_len); % Случайное сообщение

for i = 1:Data_len
if data(i) == 1
sig_inf(i*Chip_len-Chip_len+1:i*Chip_len) = chip_one;
end

if data(i) == 0
sig_inf(i*Chip_len-Chip_len+1:i*Chip_len) = chip_zero;
end

end

figure
plot(abs(sig_inf));


pream = pream/max(pream);


frame_tx = [pream, sig_inf];

figure
plot(abs(frame_tx));
Tx_Data = frame_tx.';

%% Передатчик


tx = sdrtx('Pluto','RadioID','usb:0'); % Определяем устройство передачи 
tx.ShowAdvancedProperties = true;
tx.UseCustomFilter = false;
tx.CenterFrequency = fc; % Задаем несущую частоту
tx.BasebandSampleRate = fs; % Задаем полосу частот
tx.Gain =-50; % -89.75 ... 0
transmitRepeat(tx, Tx_Data); % Передача данных Tx_Data
% tx.ShowAdvancedProperties = true;
% 

%% Приемник

rx = sdrrx('Pluto','RadioID','usb:1'); % Определяем устройство передачи
rx.CenterFrequency = fc; % Выбираем несущую частоту
rx.BasebandSampleRate = fs; % Выбираем полосу частот
rx.Gain = 10; % -4 ... 71

rx.ShowAdvancedProperties = true;
rx.EnableQuadratureCorrection = true;
rx.EnableRFDCCorrection = true;
rx.EnableBasebandDCCorrection = true;

rx.OutputDataType = 'double';
rx.SamplesPerFrame = length(Tx_Data)*Buf_len_koef; % Размер буффера



rx_data_in_time=[];
rx_data_in_time(:,1) = rx(); % Принятые данные
rx_data_in_time = rx_data_in_time.';


korr = xcorr(rx_data_in_time, pream);

[A,B] = max(korr(1:length(rx_data_in_time)+length(rx_data_in_time)/2));

start = B-length(rx_data_in_time) +1;
% sym1 = rx_data_in_time(start:start+Pream_len-1);
% sym2 = rx_data_in_time(start+Pream_len:start+2*Pream_len-1);

figure
plot(abs(korr));

sig_inf_reciv = rx_data_in_time(start+Pream_len*2:start+Pream_len*2+Chip_len*Data_len-1);


corr=xcorr(sig_inf_reciv,chip_one);
figure
plot(abs(corr));

chip_matr = reshape(sig_inf_reciv, Chip_len, Data_len);

chip_matr = chip_matr.';

for i = 1: Data_len
korr_koeff_one(i) = abs(sum(chip_matr (i,:).*conj(chip_one)));
korr_koeff_zero(i) = abs(sum(chip_matr (i,:).*conj(chip_zero)));
if (korr_koeff_one (i)>korr_koeff_zero(i))
bits_demod (i) = 1;
else
    bits_demod (i) = 0;
end
end

[number,ratio] = biterr(data,bits_demod)
