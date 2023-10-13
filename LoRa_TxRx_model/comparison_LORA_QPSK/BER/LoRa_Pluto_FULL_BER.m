clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
SF = 12;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 30e6;
fc = 2200e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;

num_pre = 8;


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
% tx_chirp = [downch, tx_preamble, mod_chirp];
tx_sync = LORA.sync;
% tx_chirp = [downch, tx_preamble, tx_sync, mod_chirp]; % DOWN, UP, SYNC, PAYLOAD
tx_chirp = [downch, tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);


%% ================================= Передача с помощью Pluto
tx = sdrtx('Pluto','RadioID','usb:0'); % 104473541196000414000800c7fa1eacf8
tx.ShowAdvancedProperties = true;
tx.UseCustomFilter = false;
tx.CenterFrequency    = fc; % Задаем несущую частоту
tx.BasebandSampleRate = BW; % Задаем полосу частот
tx.Gain = -30; %-89.75 to 0
transmitRepeat(tx, normalize(tx_chirp.'));    % Осуществляем непрерывную передачу


%% ================================= Прием с помощью Pluto
% figure(3)
for nIter=1:10

        rx = sdrrx('Pluto','RadioID','usb:1');
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
        

        %% ================================= Correlation
        [cor,lags] = xcorr(channel_chirp_sto, downch);
        cor(end/2:end) = cor(end/2:end).*gausswin(length(cor(end/2:end))).';
        [max_amp, max_idx] = max(abs(cor));
        start = lags(max_idx);
        channel_chirp_corr = channel_chirp_sto(abs(start)+1:abs(start)+tx_length);

%     for j=1:length(channel_chirp_corr)
%         channel_chirp_corr(j) = channel_chirp_corr(j).*exp(1i*DEraw*j*(-1));
%     end
%         channel_chirp_corr = channel_chirp_corr(Base+1:end);

        %% ================================= Frequency correction
%         rx_preamble = channel_chirp_corr(1:Base*num_pre);
%         corrected_signal = channel_chirp_corr(Base*num_pre+1:end);
%         [corrected_signal, rx_preamble] = LORA.CFO(channel_chirp_corr, num_pre);
        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM(channel_chirp_corr, num_pre);
%         [corrected_signal] = LORA_RLS(64, preamble, preamb_align, corrected_signal);

        %% ================================= демодуляция
        aos = 5;
        [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( corrected_signal, num_symRM, tx_preamble, rx_preamble, aos);
        
        
        %% ================================= декодирование CRC
        [data_crc_decodeRM] = LORA.decodeCRC(hard_bits, num_symRM, zeros2end, flagRM);
%         data_crc_decodeRM = hard_bits;
        
        
        %% ================================= LDPC decoding
%         maxnumiter = 100;
%         rxbits = ldpcDecode(data_crc_decodeRM.', cfgLDPCDec, maxnumiter).';
        rxbits = data_crc_decodeRM;


        %% ================================= Debugging
        bit_err = sum(abs(data~=rxbits));

        fprintf('time_shift = %d\n', start)
        fprintf('bit_err = %d\n', bit_err)
%         fprintf('est1 = %.2f\n', est1)
%         fprintf('est2 = %.2f\n', est2)
%         fprintf('est3 = %.2f\n', est3)
%         fprintf('est_full = %.2f\n', est_full)
        
        
%         stem(abs(fourier))
%         drawnow

% #########################################################################
    if(bit_err>0)
        release(tx);
        release(rx);

        % Calculation
        for mp = 2:num_pre+1
            fourier(mp-1,:) = abs(fft(channel_chirp_corr(mp*Base-Base+1:Base*mp).*downch));
            [TEraw] = LORA.TIME_ERR(fourier(mp-1,:));
            TEraw_array(mp-1) = TEraw;
        end
        TEraw_mean = mean(abs(TEraw_array));

        figure(1)
        stem(abs(fourier.'))
%         xlim([0 40])

            FEraw = -((BW*TEraw_mean)/Base)/2.685;
            DEraw = FEraw*2*pi*(1/BW);

            % Compensate in preamb
            comp_sig_1 = zeros(1,length(channel_chirp_corr));
            for j=1:length(channel_chirp_corr)
                comp_sig_1(j) = channel_chirp_corr(j).*exp(1i*DEraw*j*(-1));
                comp_sig_2(j) = channel_chirp_corr(j).*exp(1i*DEraw*j);
            end

            for mp = 2:num_pre+1
                fourier1 = abs(fft(comp_sig_1(mp*Base-Base+1:Base*mp).*downch));
                fourier2 = abs(fft(comp_sig_2(mp*Base-Base+1:Base*mp).*downch));
                [TEraw1] = LORA.TIME_ERR(fourier1);
                TEraw_array1(mp-1) = TEraw1;
                [TEraw2] = LORA.TIME_ERR(fourier2);
                TEraw_array2(mp-1) = TEraw2;
            end
            TEraw_mean1 = mean( abs(TEraw_array1) );
            TEraw_mean2 = mean( abs(TEraw_array2) );
            if(TEraw_mean1>TEraw_mean2)
                DEraw = -1*DEraw;
            end
fprintf('TEraw_mean1   = %s\n', num2str(TEraw_mean1))
fprintf('TEraw_mean2   = %s\n\n', num2str(TEraw_mean2))
% return
% % Compensation
% FEraw = -((BW*TEraw_mean)/Base)/2.685;
% DEraw = FEraw*2*pi*(1/BW);
% channel_chirp_corr = channel_chirp_realign_te;
extract_chirp = channel_chirp_corr(Base+1:Base*2);
for j=1:Base
    sig_cor(j) = extract_chirp(j).*exp(1i*DEraw*j*(-1));
end
fourier_est = abs(fft(sig_cor.*downch));
[TEraw_est] = LORA.TIME_ERR(fourier_est);
TEraw_est_reg = TEraw_est;

fprintf('TEraw_mean    = %s\n', num2str(TEraw_mean))
fprintf('TEraw_array    = %s\n', num2str(TEraw_array))
fprintf('TEraw_est_reg = %s\n', num2str(TEraw_est_reg))

for k=1:100
end
% return
figure(2)
% plot((Ri))
stem(abs(fourier_est))
% xlim([0 50])
return
        err_pos = find(check_data_no_gray~=sv_cor);
        fprintf('check  = %s\n', num2str(check_data_no_gray(err_pos)))
        fprintf('sv_cor = %s\n', num2str(sv_cor(err_pos)))

                break
% %% ================================= 3. Golden rectangle опт-ая версия
% accu = 1;                                  % Точность оценки частоты
% r = channel_chirp_realign(1:Base*num_pre);  % Входной сигнал (один чирп)
% C_conj = repmat(downch,1,num_pre);          %
% fdd = ceil(fps/2);                          % Граница диапазона ошибок
% RD = fdd*2;                                 % Диапазон ошибок
% k=1:Base;
% 
% % Золотое сечение
% for j=1:accu:RD
%     
%      % защита от переполнения
%     if (j>(RD))
%     j=RD;
%     end
% 
% % Первое приближение
% df(j) = -fdd+j;
% 
% % формируем 8 опорных чирпов с компенсацией
% for n=1:num_pre
% Local_chirp(n*Base-Base+1:Base*n) = exp(-1i*2*pi*(df(j)*(k+n*Base))/BW);
% end
% 
% % Компенсируем у преамбулы сдвиг и находим корр функцию
% Zi = r.*Local_chirp;
% Ri(j) =  abs(sum((Zi.*C_conj)))/num_pre;
% end

    end
% #########################################################################

    % Estimate the BER for both methods
    BER(nIter) = bit_err;
%     figure(2)
%     subplot(211)
%     plot(lags, abs(cor))
%     subplot(212)
%     plot(lags)
end

release(tx);
release(rx);
fprintf('num_err = %s\n', num2str(BER))
fprintf('mean_err = %s\n', num2str(mean(BER)))



% save('rx_data_in_time.mat','rx_data_in_time')
% save('data.mat','data')

% function [TEraw] = TIME_ERR(product, Base);
%     
%     [~, ind_max] = max( abs(product) );
%     if((ind_max-1)<=0)
%         mag1 = abs(product(Base+(ind_max-1)));
%     else
%         mag1 = abs(product(ind_max-1));
%     end
% 
%     mag2 = abs(product(ind_max));
% 
%     if((ind_max+1)>Base)
%         mag3 = abs(product((ind_max+1)-Base));
%     else
%         mag3 = abs(product(ind_max+1));
%     end
%     TEraw = (mag3+mag1)/mag2;
% 
% %     if((ind_max+1)>=0)
% %         mag3 = abs(product(ind_max+1));
% %     else
% %         mag3 = abs(product(ind_max+1));
% %     end
% %     TEraw = (mag3-mag1)/mag2;
% end