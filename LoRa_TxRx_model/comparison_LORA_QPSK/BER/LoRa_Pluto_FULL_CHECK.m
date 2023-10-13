clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
SF = 11;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 30e6;
fc = 920e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;

num_pre = 8;
warning('off')

%% ================================= LDPC coding
% [cfgLDPCEnc,cfgLDPCDec] = LORA.generateConfigLDPC(3/4, 1944);
%         
% % Number of message bits
% numinfobits = cfgLDPCEnc.NumInformationBits; 
% numcodebits = cfgLDPCEnc.BlockLength; 
% 
% % Message/Iformation bits
% data = randi([0 1],1, numinfobits); 
% data_ldpc_code = ldpcEncode(data.', cfgLDPCEnc).';

% numcodebits = 1944;
numcodebits = 648;
data = randi([0 1],1, numcodebits); 
data_ldpc_code = data;

%% ================================= Rate matching
[data_ldpc_codeRM, numcodebitsRM, num_symRM, zeros2end, flagRM] = LORA.RM(data_ldpc_code);


%% ================================= CRC coding
[data_crc_ldpc_codeRM] = LORA.codeCRC(data_ldpc_codeRM, num_symRM);
% data_crc_ldpc_codeRM = data_ldpc_codeRM;
% num_symRM = numcodebits/SF;


%% ================================= Mодуляция
[mod_chirp, check_data, check_data_no_gray] = LORA.lorax_modified( data_crc_ldpc_codeRM, num_symRM, 1);
tx_preamble = repmat(LORA.chirp,1,num_pre);
tx_chirp = [downch, tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);
sig_length = tx_length-Base;


%% ================================= Передача с помощью Pluto
% tx = sdrtx('Pluto','RadioID','usb:0'); % 104473541196000414000800c7fa1eacf8
load mycustomfilter
% filtnv{2} = 30e6;
% filtnv{14} = 30e6;
% filtnv_{8} = 16;
% filtnv_{6} = 10000*ones(1,16);
% filtnv{10} = 4;
% filtnv{12} = 0;
tx = sdrtx('Pluto','RadioID','usb:0',filtnv{:});
tx.UseCustomFilter = true;

tx.ShowAdvancedProperties = true;
tx.UseCustomFilter = false;
tx.CenterFrequency    = fc; % Задаем несущую частоту
tx.BasebandSampleRate = BW; % Задаем полосу частот
tx.Gain = -60; %-89.75 to 0
transmitRepeat(tx, normalize(tx_chirp.'));    % Осуществляем непрерывную передачу


%% ================================= Прием с помощью Pluto
figure(3)
for nIter=1:10

%         rx = sdrrx('Pluto','RadioID','usb:0');
        load mycustomfilter_rx
%         filtnv_rx{2} = 30e6;
%         filtnv_rx{14} = 30e6;
%         filtnv_rx{8} = 16;
%         filtnv_rx{6} = 10000*ones(1,16);
%         filtnv_rx{10} = 4;
%         filtnv_rx{12} = 0;
        rx = sdrrx('Pluto','RadioID','usb:1',filtnv_rx{:});
        rx.UseCustomFilter = true;

        % info(rx) %104473541196000508002e000e2d1ffaa3
        rx.CenterFrequency    = fc; % Выбираем несущую частоту
        rx.BasebandSampleRate = BW; % Выбираем полосу частот
        rx.SamplesPerFrame = tx_length*5; % Принимаем в 10 раз больше чем отправили
        
        % Включаем слежение за квадратурами
        rx.ShowAdvancedProperties = true;
        rx.EnableQuadratureCorrection = true;
        rx.EnableRFDCCorrection = true;
        rx.EnableBasebandDCCorrection = true;
        rx.OutputDataType = 'double';
        rx.GainSource = 'AGC Slow Attack';
        
        rx_data_in_time = rx();
        channel_chirp_sto = rx_data_in_time.';
%         channel_chirp_sto = resample(channel_chirp_sto, 1, rsfac);

        %% ================================= Correlation
        [cor,lags] = xcorr(channel_chirp_sto, downch);
        cor(round(end/2):end) = cor(round(end/2):end).*gausswin(length(cor(round(end/2):end))).';
        [max_amp, max_idx] = max(abs(cor));
        start = lags(max_idx);
%         channel_chirp_corr = channel_chirp_sto(abs(start)+1:abs(start)+tx_length);
        channel_chirp_corr = channel_chirp_sto(abs(start)+1:end);

%     for j=1:length(channel_chirp_corr)
%         channel_chirp_corr(j) = channel_chirp_corr(j).*exp(1i*DEraw*j*(-1));
%     end

%         channel_chirp_corr = filtfilt(ones(1, 8)/(8), 1, channel_chirp_corr);
%         channel_chirp_corr = filtfilt(d1,channel_chirp_corr);


        %% ================================= Frequency correction
%         channel_chirp_corr = channel_chirp_corr(Base+1:end);
%         rx_preamble = channel_chirp_corr(1:Base*num_pre);
%         corrected_signal = channel_chirp_corr(Base*num_pre+1:end);
%         [freq_data, corrected_signal, rx_preamble] = LORA.CFO(channel_chirp_corr, num_pre);
        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM(channel_chirp_corr, num_pre, sig_length);
%         [corrected_signal] = LORA_RLS(14, tx_preamble, rx_preamble, corrected_signal);

% figure(1)
% plot(real(tx_preamble))
% hold on
% plot(real(rx_preamble))
% release(tx);
% release(rx);
% return
        %% ================================= демодуляция
        aos = 2;
        [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( corrected_signal, num_symRM, tx_preamble, rx_preamble, aos);
%         [sv_decode, fourier] = LORA.delorax_cyclic(num_symRM, corrected_signal);
%         dbits = de2bi(sv_decode, SF)';
%         hard_bits = dbits(:)';
        
        %% ================================= декодирование CRC
        [data_crc_decodeRM] = LORA.decodeCRC(hard_bits, num_symRM, zeros2end, flagRM);
%         data_crc_decodeRM = hard_bits;
        
        
        %% ================================= LDPC decoding
%         maxnumiter = 100;
%         rxbits = ldpcDecode(data_crc_decodeRM.', cfgLDPCDec, maxnumiter).';
        rxbits = data_crc_decodeRM;


        %% ================================= Debugging
        bit_err = sum(abs(data~=rxbits));

        STOint = freq_data{1};
        est1 = freq_data{2};
        est4 = freq_data{3};
        FEraw = freq_data{4};
        
        
        fprintf('time_shift = %d\n', start)
        fprintf('bit_err = %d\n', bit_err)
        fprintf('STOint = %.6f\n', STOint)
        fprintf('est1 = %.6f\n', est1)
        fprintf('est4 = %.6f\n', est4)
        fprintf('FEraw = %.6f\n', FEraw)

        stem(abs(fourier))
        drawnow

%         if(bit_err>0)
%             break
%         end
    % Estimate the BER for both methods
    BER(nIter) = bit_err;
end

release(tx);
release(rx);



% save('rx_data_in_time.mat','channel_chirp_sto')
% save('data.mat','data')

fprintf('num_err = %s\n', num2str(BER))
fprintf('mean_err = %s\n', num2str(mean(BER)))


return
figure(10)
subplot(221)
stem(abs(fft(rx_preamble(Base+1:2*Base).*downch)))
title('before correction')
subplot(222)
plot(real(rx_preamble(Base+1:2*Base)))
title('before correction')
subplot(223)
stem(abs(fft(corrected_signal(1:Base).*downch)))
title('after correction')
subplot(224)
plot(real(corrected_signal(1:Base)))
title('after correction')

% figure(1)
% plot(sv-check_data)
function [TEraw] = TIME_ERR(product, Base);
    
    [~, ind_max] = max( abs(product) );
    if((ind_max-1)<=0)
        mag1 = abs(product(Base+(ind_max-1)));
    else
        mag1 = abs(product(ind_max-1));
    end

    mag2 = abs(product(ind_max));

%     if((ind_max+1)>=0)
%         mag3 = abs(product(ind_max+1));
%     else
%         mag3 = abs(product(ind_max+1));
%     end
    TEraw = (mag3-mag1)/mag2;
end