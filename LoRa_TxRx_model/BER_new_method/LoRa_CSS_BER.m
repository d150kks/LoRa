clc
clear all
close all

%% ================================= Переменные

% коэффициенты
SF = 8;        % коэффициент расширения спектра (от 7 до 12)
Base = 2^SF;   % База сигнала
num_sym = 1;% Число передаваемых символов
num_pre = 8;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
snr = [-20:2:0];

% массивы данных
numbits = SF*num_sym;          % число бит
% data = randi([0, 1], 1, numbits); % формирование массива бит
data= [1	1	0	0	0	0	0	0 1	1	0	0	0	0	0	0];

% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 125e3+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 1000e6;            % центральная частота
band = BW;
Ts = (2^SF)/BW;        % длительность сигнала


%% ================================= модуляция
os = 1;
fs = BW*os;
ts=1/fs;

data_rep = zeros(1,numbits*2);
for i=1:num_sym
    data_rep(SF*i-SF+1:SF*i) = data(SF*i-SF+1:SF*i);
    data_rep(SF*(i+1)-SF+1:SF*(i+1)) = data(SF*i-SF+1:SF*i);
end

[mod_chirp, check_chirp, downch, check_data] = lorax_modified( data_rep, fs, Base, Ts, SF, BW);
N = length(downch);


tx_length = length(mod_chirp);
% Forming Signal
freq_shift = 800;

% вводим частотный сдвиг
dphi=freq_shift*2*pi*ts;% сдвиг

% add freq shift
for j=1:tx_length
  tx_chirp(j)=mod_chirp(j)*exp(1i*dphi*j);
end





%% ================================= BER
tic
for n = 1:length(snr)

    [numErr, NumData] = deal(0);

    while  NumData < 1e5
        % АБГШ демодуляция и декодирование
        snr(n)=-14;
        rxSig = awgn(tx_chirp,snr(n),'measured');
        num = length(downch);
B=2^SF;
% css = [];
i=1;

d = rxSig(num*i-num+1:num*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
d2 = rxSig(num*(i+1)-num+1:num*(i+1)).*downch;

fourier = abs(fft(d));            % переводим результат в область частот
fourier2 = abs(fft(d2));
% fourier = abs(fourier(end/2:end));
[~, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
[~, indexMax2] = max( fourier2 ); % находим щелчок  частоты в чирпе
% вычисляем значение кодового слова исходя из базы сигнала

%         [sv, fourier] = delorax_modified( length(data), SF, downch, rxSig);
%         dbits = de2bi(sv,SF)';
%         dbits = dbits(:)';
        
figure(1)
subplot(211)
stem(abs(fourier))
hold on 
stem(abs(fourier2))
xlim([0, 40])
subplot(212)
stem(abs(fourier+(fourier2)))
title('sum')
xlim([0, 40])

ff = fft(d2);
angle(ff(4))/(2*pi*Ts)

% figure(1)
% stem(angle(ff))
% xlim([0, 40])


        return
        % подсчет БЕР с учетом задержки
        err = sum(dbits~=data);

        % Increment the error and bit counters
        numErr = numErr + err;        
        NumData = NumData + numbits;
    end

    % Estimate the BER for both methods
    BER(n) = numErr/NumData;
end
toc

%%
figure(1)
semilogy(snr,BER,'-*','color','k');
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');

% save('lora_css_ber.mat','BER')
% save('lora_css_snr.mat','snr')


function [sv, fourier] = delorax_modified2( length_data, SF, downch, chirp)
num = length(downch);
B=2^SF;
% css = [];
for i = 1:(length_data/SF)-1

d = chirp(num*i-num+1:num*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
d2 = chirp(num*(i+1)-num+1:num*(i+1)).*downch;

fourier = abs(fft(d));            % переводим результат в область частот
fourier2 = abs(fft(d2));
% fourier = abs(fourier(end/2:end));
[~, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
[~, indexMax2] = max( fourier2 ); % находим щелчок  частоты в чирпе
% вычисляем значение кодового слова исходя из базы сигнала

% sv(i) = indexMax-1;
if indexMax>B
    sv(i) = B - (num - indexMax) - 1;
else
    sv(i) = indexMax - 1;
end

end

end