clc
clear all
close all


tic

%% ================================= Переменные
fc = 900e6;
fs = 30e6;

M = 4;   
Pream_len = 4096;
Chip_len = 300;
Data_len = 648;
warning('off')

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


%% Forming Signal
tx_sig = [pream, sig_mod];
tx_length = length(tx_sig);

%% ================================= Передача с помощью Pluto
% tx = sdrtx('Pluto','RadioID','usb:0'); % 104473541196000414000800c7fa1eacf8
FAMP = 1100;
load mycustomfilter
tx = sdrtx('Pluto','RadioID','usb:0',filtnv{:});
tx.UseCustomFilter = true;
tx.filtCoefficients = [FAMP*32, zeros(1,31)];
tx.filtGain = 0;

tx.ShowAdvancedProperties = true;
tx.CenterFrequency    = fc; % Задаем несущую частоту
tx.BasebandSampleRate = fs; % Задаем полосу частот
tx.Gain = -20; %-89.75 to 0
transmitRepeat(tx, normalize(tx_sig.'));    % Осуществляем непрерывную передачу
        
% release(tx)
% return
%% ================================= Прием с помощью Pluto
% figure(3)
for nIter=1:10

%         rx = sdrrx('Pluto','RadioID','usb:1');
        load mycustomfilter_rx
        rx = sdrrx('Pluto','RadioID','usb:1',filtnv_rx{:});
        rx.UseCustomFilter = true;
        rx.filtCoefficients = [FAMP*32, zeros(1,31)];
        rx.filtGain = 0;

        rx.CenterFrequency    = fc; % Выбираем несущую частоту
        rx.BasebandSampleRate = fs; % Выбираем полосу частот
        rx.SamplesPerFrame = tx_length*5; % Принимаем в 10 раз больше чем отправили
        
        % Включаем слежение за квадратурами
        rx.ShowAdvancedProperties = true;
        rx.EnableQuadratureCorrection = true;
        rx.EnableRFDCCorrection = true;
        rx.EnableBasebandDCCorrection = true;
        rx.OutputDataType = 'double';
%         rx.GainSource = 'AGC Slow Attack';
        rx.GainSource = 'Manual';
        rx.Gain = 50;

        rx_data_in_time = rx();
        rx_sig = rx_data_in_time.';
        

        %% ================================= Correlation
        [cor,lags] = xcorr(rx_sig, pream);
        cor(end/2:end) = cor(end/2:end).*gausswin(length(cor(end/2:end))).';
        [max_amp, max_idx] = max(abs(cor));
        start = lags(max_idx);
        rx_sig_corr = rx_sig(abs(start)+Pream_len*2+1:abs(start)+tx_length);

% figure(1)
% plot(real(sig_mod))
% hold on
% plot(real(rx_sig_corr))
% % xlim([1 Base*5])
% return
        %% ================================= Demodulation
        demod_data = QPSK_DEMOD(rx_sig_corr, chip_one, chip_zero, Chip_len, Data_len);

        %% ================================= Debugging
        bit_err = sum(abs(data-demod_data));

        % Increment the error and bit counters

        fprintf('time_shift = %d\n', start)
        fprintf('bit_err = %d\n', bit_err)

        
%         stem(abs(fourier))
%         drawnow

    % Estimate the BER for both methods
    BER(nIter) = bit_err;

end

release(tx);
release(rx);
fprintf('num_err = %s\n', num2str(BER))
fprintf('mean_err = %s\n', num2str(mean(BER)))

% figure(3)
% stem(abs(fourier))
% 
% figure(1)
% semilogy(tx_power,BER,'-*','color','k');
% grid
% xlabel('SNR (dB)')
% ylabel('Bit Error Rate')
% hold off
% title('SNR');

% 
% save('lora_crc_pluto_ber.mat','BER')
% save('lora_crc_pluto_snr.mat','tx_power')

%% ================================= Demodulation
function demod_data = QPSK_DEMOD(rxSig, chip_one, chip_zero, Chip_len, Data_len)

    
    win_one = zeros(1,Chip_len);
    win_zero = zeros(1,Chip_len);
    demod_data = zeros(1,Data_len);
    for i=1:Data_len
        win_one = abs(sum(rxSig(i*Chip_len-Chip_len+1:i*Chip_len).*conj(chip_one)));
        win_zero = abs(sum(rxSig(i*Chip_len-Chip_len+1:i*Chip_len).*conj(chip_zero)));
        if(win_one>win_zero)
            demod_data(i) = 1;
        else
            demod_data(i) = 0;
        end
    end

end