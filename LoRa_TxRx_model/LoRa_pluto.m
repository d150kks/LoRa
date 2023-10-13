clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
SF = 8;        % коэффициент расширения спектра (от 7 до 12)
BW = 30e6;
LORA = myLoRaClass(SF, BW);
chirp = LORA.chirp;
Base = LORA.Base;
downch = LORA.downch;

num_sym = 100;% Число передаваемых символов
num_pre = 8;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
% snr = 40;      % ОСШ

% массивы данных
numbits = SF*num_sym;          % число бит
data = randi([0, 1], 1, numbits); % формирование массива бит

% частота, фаза, время
fc = 2200e6;            % центральная частота
Ts = (2^SF)/BW;        % длительность сигнала


%% ================================= модуляция
os = 1;
fs = BW*os;
ts=1/fs;

% ##################
data_resh = reshape(data, SF, []);
data_rep = repmat(data_resh, 1, 1);
data_rep = reshape(data_rep, 1, []);
% ##################

% [mod_chirp, check_chirp, downch, check_data] = lorax_modified( data_rep, fs, Base, Ts, SF, BW);
[mod_chirp, check_data] = LORA.lorax_modified( data, num_sym, 1);
N = length(downch);

%% Forming Signal
preamble = repmat(chirp,1,num_pre);
%%tx_chirp = [preamble, mod_chirp];
tx_chirp = [downch, preamble, mod_chirp];

tx_length = length(tx_chirp);
sig_length = tx_length-Base;



% return
%% ================================= Передача с помощью Pluto
FAMP = 1100;
load mycustomfilter
tx = sdrtx('Pluto','RadioID','usb:0',filtnv{:});
tx.UseCustomFilter = true;
tx.filtCoefficients = [FAMP*32, zeros(1,31)];
tx.filtGain = 0;

tx.ShowAdvancedProperties = true;
tx.CenterFrequency    = fc; % Задаем несущую частоту
tx.BasebandSampleRate = BW; % Задаем полосу частот
tx.Gain = -35; %-89.75 to 0
transmitRepeat(tx, normalize(tx_chirp.'));    % Осуществляем непрерывную передачу

%% ================================= Прием с помощью Pluto
figure(3)
for nIter=1:10

        load mycustomfilter_rx
        rx = sdrrx('Pluto','RadioID','usb:1',filtnv_rx{:});
        rx.UseCustomFilter = true;
        rx.filtCoefficients = [FAMP*32, zeros(1,31)];
        rx.filtGain = 0;

        rx.CenterFrequency    = fc; % Выбираем несущую частоту
        rx.BasebandSampleRate = BW; % Выбираем полосу частот
        rx.SamplesPerFrame = tx_length*10; % Принимаем в 10 раз больше чем отправили
        
        % Включаем слежение за квадратурами
        rx.ShowAdvancedProperties = true;
        rx.EnableQuadratureCorrection = true;
        rx.EnableRFDCCorrection = true;
        rx.EnableBasebandDCCorrection = true;
        rx.OutputDataType = 'double';
        rx.GainSource = 'Manual';
        rx.Gain = 50;
        
        rx_data_in_time = rx();
        channel_chirp_sto = rx_data_in_time.';
        
        
        %% ================================= Correlation
        [cor,lags] = xcorr(channel_chirp_sto, downch);
        cor(round(end/2):end) = cor(round(end/2):end).*gausswin(length(cor(round(end/2):end))).';
        [max_amp, max_idx] = max(abs(cor));
        start = lags(max_idx);
        channel_chirp_corr = channel_chirp_sto(abs(start)+1:end);

        %% ================================= Frequency correction
%         [freq_data, corrected_signal] = LORA.LORA_FREQ_ESTIM(channel_chirp_corr, num_pre);
        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM(channel_chirp_corr, num_pre, sig_length);
%         [corrected_signal] = LORA_RLS(64, preamble, preamb_align, corrected_signal);

%         figure(1)
%         plot(real(mod_chirp))
%         hold on
%         plot(real((corrected_signal)))
%         return
        %% ================================= Demodulation
        [sv_decode, sv, fourier] = LORA.delorax_modified(corrected_signal, num_sym);
%         [sv, fourier] = delorax_modified( length(data), SF, downch, corrected_signal);
%         [sv, fourier] = LORA.delorax_cyclic(num_sym, corrected_signal);
        dbits = de2bi(sv_decode,SF).';
        dbits = dbits(:).';


% num = length(downch);
% B=2^SF;
% % css = [];
% for i = 1:length(data_rep)/SF
% 
% pair_chirp = corrected_signal((i-1)*num*3+1: num*3*(i));
% d = pair_chirp(1:num).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
% d2 = pair_chirp(num+1:num*2).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
% d3 = pair_chirp(num*2+1:num*3).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
% 
% fourier = abs(fft(d));            % переводим результат в область частот
% fourier2 = abs(fft(d2));            % переводим результат в область частот
% fourier3 = abs(fft(d3));            % переводим результат в область частот
% fourier4 = abs(fft( [pair_chirp(1:num), pair_chirp(num+1:num*2), pair_chirp(num*2+1:num*3)].*[downch downch downch]  ));
% [~, indexMax] = max( fourier+fourier2 ); % находим щелчок  частоты в чирпе
% % css = fourier
% %%%%%%%
% % вычисляем значение кодового слова исходя из базы сигнала
% 
% % sv(i) = indexMax-1;
% if indexMax>B
%     sv(i) = B - (num - indexMax) - 1;
% else
%     sv(i) = indexMax - 1;
% end
% break
% end
% 
% figure(10)
% subplot(211)
% stem(abs(fourier))
% hold on
% stem(abs(fourier2))
% 
% subplot(212)
% stem( fourier+fourier2+fourier3 )
% 
% figure(11)
% stem(fourier4)
% return

        %% ================================= Debugging
        bit_err = sum(abs(data-dbits));

        est1 = freq_data{1};
        est2 = freq_data{2};
        est3 = freq_data{3};
        est_full = freq_data{4};
        
        fprintf('time_shift = %d\n', start)
        fprintf('bit_err = %d\n', bit_err)
        fprintf('est1 = %.2f\n', est1)
        fprintf('est2 = %.2f\n', est2)
        fprintf('est3 = %.2f\n', est3)
        fprintf('est_full = %.2f\n', est_full)
        
        
        stem(abs(fourier))
        drawnow

    % Estimate the BER for both methods
    BER(nIter) = bit_err;

end

release(tx);
release(rx);
fprintf('num_err = %s\n', num2str(BER))
fprintf('mean_err = %s\n', num2str(mean(BER)))



% figure(1)
% plot(data_code_polar-data_decode, 'LineWidth',2)
% hold on
% % plot(data_decode)
% plot(data_decode_polar-data, 'LineWidth',2)
% legend('err before polar decoding', 'err after polar decoding')
% % plot(data_decode_polar, 'LineWidth',2)
% % hold on
% % plot(data)



function [sv, fourier] = delorax_modified2( length_data, SF, downch, chirp)
num = length(downch);
B=2^SF;
% css = [];
for i = 1:(length_data/SF)/1

pair_chirp = chirp((i-1)*num*4+1: num*4*(i));
d = pair_chirp(1:num).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
d2 = pair_chirp(num+1:num*2).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
d3 = pair_chirp(num*2+1:num*3).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
d4 = pair_chirp(num*3+1:num*4).*downch;   % перемножаем входной и опорный ОБРАТНый чирп

fourier = (fft(d));            % переводим результат в область частот
fourier2 = (fft(d2));            % переводим результат в область частот
fourier3 = (fft(d3));            % переводим результат в область частот
fourier4 = (fft(d4));            % переводим результат в область частот
[~, indexMax] = max( abs(fourier+fourier2+fourier3+fourier4) ); % находим щелчок  частоты в чирпе
% css = fourier
%%%%%%%
% вычисляем значение кодового слова исходя из базы сигнала

sv(i) = indexMax-1;
% indexMax-1;
% if indexMax>B
%     sv(i) = B - (num - indexMax) - 1;
% else
%     sv(i) = indexMax - 1;
% end

end

end


