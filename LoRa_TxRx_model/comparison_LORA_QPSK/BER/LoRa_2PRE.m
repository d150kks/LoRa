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

num_pre = 2;


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

%% ================================= Signal Constructing
LORA2 = myLoRaClass(SF+2,BW);
Base2 = LORA2.Base;
downch2 = LORA2.downch;
upch2 = LORA2.chirp;
tx_preamble = [upch2, upch2];
tx_chirp = [downch2, tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);

% return
%% ================================= Передача с помощью Pluto
tx = sdrtx('Pluto','RadioID','usb:0'); % 104473541196000414000800c7fa1eacf8
tx.ShowAdvancedProperties = true;
tx.UseCustomFilter = false;
tx.CenterFrequency    = fc; % Задаем несущую частоту
tx.BasebandSampleRate = BW; % Задаем полосу частот
tx.Gain = -40; %-89.75 to 0
transmitRepeat(tx, normalize(tx_chirp.'));    % Осуществляем непрерывную передачу


%% ================================= Прием с помощью Pluto
figure(3)
for nIter=1:100

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
%         channel_chirp_sto = resample(channel_chirp_sto, 1, rsfac);

        %% ================================= Correlation
        [cor,lags] = xcorr(channel_chirp_sto, downch2);
        cor(end/2:end) = cor(end/2:end).*gausswin(length(cor(end/2:end))).';
        [max_amp, max_idx] = max(abs(cor));
        start = lags(max_idx);
%         channel_chirp_corr = channel_chirp_sto(abs(start)+1:abs(start)+tx_length);
        channel_chirp_corr = channel_chirp_sto(abs(start)+1:end);


        %% ================================= Frequency correction
%         N = Base2;
%         fps = BW/Base2; % Определяем Hz/samp
%         rx_downch = channel_chirp_corr(1:N);
%         rx_preamb = channel_chirp_corr(N+1:N*3);
% 
%         % STO and CFO estimation
%         [~, ~, fup_idx, fdown_idx, ~, ~] = LORA2.fraq(rx_preamb, rx_downch);
%         STOint = ceil((fup_idx-fdown_idx)/2);
% 
%         channel_chirp_corr = channel_chirp_corr((1+N)-STOint:end-STOint);
%         channel_chirp_treshold = channel_chirp_corr(1:num_pre*N);
% 
%         % ~~~~~~~~ 1. Coarse estimation ~~~~~~~~ 
%         for i = 1:num_pre
%             fourier = abs(fft(channel_chirp_treshold(i*N-N+1:N*i).*downch2));
%             [~, ind1] = max( fourier );
%             pre_align(i) = ind1;
%             
% %             if (pre_align(i)>N/2)
% %                 pre_align(i) = pre_align(i)-N;
% %             end
%         end
% % 
% %         figure(1)
% %         stem(fourier)
% %         return
%         est1 = (mean(pre_align-1))*fps;
%         dphi1 = est1*2*pi*(1/BW); % сдвиг
% 
%         for j = 1:length(channel_chirp_corr)
%             corrected_signal(j) = channel_chirp_corr(j).*exp(1i*dphi1*j*(-1));
%         end
%         rx_preamble = corrected_signal(1:num_pre*N);
%         corrected_signal = corrected_signal(num_pre*N+1:end);
        num_pre = 2;
        sig_length = tx_length-length(downch2);
        [freq_data, corrected_signal, rx_preamble] = LORA2.LORA_FREQ_ESTIM(channel_chirp_corr, num_pre, sig_length);
%         [corrected_signal] = LORA_RLS(24, tx_preamble, rx_preamble, corrected_signal);

        %% ================================= демодуляция
%         aos = 5;
%         [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( corrected_signal, num_symRM, tx_preamble, rx_preamble, aos);
        [sv_decode, fourier] = LORA.delorax_cyclic(num_symRM, corrected_signal);
        dbits = de2bi(sv_decode, SF)';
        hard_bits = dbits(:)';
        
        %% ================================= декодирование CRC
        [data_crc_decodeRM] = LORA.decodeCRC(hard_bits, num_symRM, zeros2end, flagRM);
%         data_crc_decodeRM = hard_bits;
        
        
        %% ================================= LDPC decoding
%         maxnumiter = 100;
%         rxbits = ldpcDecode(data_crc_decodeRM.', cfgLDPCDec, maxnumiter).';
        rxbits = data_crc_decodeRM;


        %% ================================= Debugging
        bit_err = sum(abs(data~=rxbits));

%         est1 = freq_data{1};
%         est2 = freq_data{2};
%         est3 = freq_data{3};
%         est_full = freq_data{4};
        
        
        fprintf('time_shift = %d\n', start)
        fprintf('bit_err = %d\n', bit_err)
%         fprintf('est1 = %.2f\n', est1)
%         fprintf('est2 = %.2f\n', est2)
%         fprintf('est3 = %.2f\n', est3)
%         fprintf('est_full = %.2f\n', est_full)

        stem(abs(fourier))
        drawnow

        if(bit_err>0)
            break
        end
    % Estimate the BER for both methods
    BER(nIter) = bit_err;
end

release(tx);
release(rx);
fprintf('num_err = %s\n', num2str(BER))
fprintf('mean_err = %s\n', num2str(mean(BER)))



save('rx_data_in_time.mat','channel_chirp_sto')
save('data.mat','data')


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