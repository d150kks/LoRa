clc
clear all
close all

%% ================================= Переменные

% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
Base = 2^SF;   % База сигнала
num_sym = 1;% Число передаваемых символов
num_pre = 8;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
% snr = [-20:2:0];

% массивы данных
numbits = SF*num_sym;          % число бит
% data = randi([0, 1], 1, numbits); % формирование массива бит
data= [1	1	0	0	0	0	0	0];

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

[mod_chirp, check_chirp, downch, check_data] = lorax_modified( data, fs, Base, Ts, SF, BW);
N = length(downch);



%% Form preamble
% sig = [check_chirp, check_chirp];
sig = [downch, mod_chirp, check_chirp];
tx_length = length(sig);


%% Add freq offset and AWGN
freq_shift = 600;
dphi=freq_shift*2*pi*ts;% сдвиг

% add freq shift
for j=1:tx_length
  tx_chirp(j)=sig(j)*exp(1i*dphi*j);
end

snr = -5;
delay = 60;
tx_chirp = awgn([zeros(1,delay) tx_chirp zeros(1,1000)], snr, 'measured');
% tx_chirp = [0 0 0 tx_chirp(1:end-3)];


fourier1 = fft(tx_chirp(N+1:2*N).*(downch));
[MaxVal1, indexMax1] = max( abs(fourier1) );
TE = (abs(fourier1(indexMax1-1))-abs(fourier1(indexMax1+1)))/abs(fourier1(indexMax1))
% fourier1(abs(fourier1)==MaxVal1) = 0;
% [MaxVal2, indexMax2] = max( abs(fourier1) );

% mod(MaxVal1-MaxVal2,N)

figure(2)
stem(abs(fourier1))
% hold on
% stem(abs(fourier2))
return
%% ================================= Correlation
[cor,lags] = xcorr(downch, tx_chirp);
% cor(1:end/2) = cor(1:end/2).*gausswin(length(cor(1:end/2))).';
[max_amp, max_idx] = max(abs(cor));
start = lags(max_idx);
% if(start<0)
  channel_chirp_corr = tx_chirp(abs(start)+1:abs(start)+tx_length);
% elseif(start>0)
%   channel_chirp_corr = tx_chirp(start:tx_length-start);
% end

% figure(1)
% subplot(211)
% plot(lags, abs(cor))
% subplot(212)
% plot(lags)
% return


%% STO Estimation
fourier2 = fft(channel_chirp_corr(1:N).*conj(downch));
fourier1 = fft(channel_chirp_corr(N+1:2*N).*(downch));

[~, indexMax1] = max( abs(fourier1) );
fup = (indexMax1-1);
phi1 = angle(fourier1(indexMax1));
[~, indexMax2] = max( abs(fourier2) );
fdown = (indexMax2-1);
phi2 = angle(fourier2(indexMax2));

% figure(2)
% stem(abs(fourier1))
% hold on
% stem(abs(fourier2))
% xlim([0, 30])

% STOint = N/2-(fup-fdown)*0.5-1
STOint = (fup-fdown)*0.5;

channel_chirp_corr2 = [channel_chirp_corr(STOint:end), zeros(1,STOint)];

% CFOfrac = BW*(phi2-phi1)/(2*pi*N)
%% CFO Estimation
fprintf('delay = %d\n', delay)
fprintf('start = %d\n', start)
fprintf('STOint = %d\n', STOint)
% fprintf('freq_shift = %.2f\n', freq_shift)
% fprintf('freq_err = %.2f\n', freq_shift-1.5*975+CFOfrac)

% fprintf('bit_err = %d\n', bit_err)
% fprintf('est1 = %.2f\n', est1)
% fprintf('est2 = %.2f\n', est2)
% fprintf('est3 = %.2f\n', est3)
% fprintf('est_full = %.2f\n', est_full)
return

% fourier1 = fft(tx_chirp(1:N).*conj(B));
% fourier2 = fft(tx_chirp(N+1:2*N).*(B));


% CFOint = ((fup+fdown)/2)
% 


return

% fourier = [];
% fourier2 = [];
% for i=1:10
% %     fourier = [fourier, fft(rx_preamb(i*N-N+1:i*N))];
% %     fourier2 = [fourier2, fft(rx_preamb2(i*N-N+1:i*N))];
%     fourier = [fourier, fft(two.*twoup)];
%     fourier2 = [fourier2, fft(two)];
% end

% fourier = fft(two.*twoup);

figure(1)
plot(real(up7))
plot(imag(up7))
hold on
% plot(abs(fourier2))
% tx_length = length(mod_chirp);
% % Forming Signal
% freq_shift = 0;
% 
% % вводим частотный сдвиг
% dphi=freq_shift*2*pi*ts;% сдвиг
% 
% % add freq shift
% for j=1:tx_length
%   tx_chirp(j)=mod_chirp(j)*exp(1i*dphi*j);
% end





%% ================================= BER
% snr=-10;
% rxSig = awgn(tx_chirp,snr,'measured');



return


figure(2)
subplot(211)
stem(abs(fourier))
hold on
stem(abs(fourier2))
xlim([0, 40])
subplot(212)
stem(abs(fourier+conj(fourier2)))
xlim([0, 40])


ff = fft(d2);
angle(ff(4))/(2*pi*Ts)

figure(1)
stem(angle(ff))
xlim([0, 40])


        return
        % подсчет БЕР с учетом задержки
        err = sum(dbits~=data);



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


% function [sv, fourier] = delorax_modified2( length_data, SF, downch, chirp)
% num = length(downch);
% B=2^SF;
% % css = [];
% for i = 1:(length_data/SF)-1
% 
% d = chirp(num*i-num+1:num*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
% d2 = chirp(num*(i+1)-num+1:num*(i+1)).*downch;
% 
% fourier = abs(fft(d));            % переводим результат в область частот
% fourier2 = abs(fft(d2));
% % fourier = abs(fourier(end/2:end));
% [~, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
% [~, indexMax2] = max( fourier2 ); % находим щелчок  частоты в чирпе
% % вычисляем значение кодового слова исходя из базы сигнала
% 
% % sv(i) = indexMax-1;
% if indexMax>B
%     sv(i) = B - (num - indexMax) - 1;
% else
%     sv(i) = indexMax - 1;
% end
% 
% end
% 
% end