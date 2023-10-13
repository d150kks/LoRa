clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
SF = 11;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 30e6;
fc = 900e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;

intrlv_state = 13;
num_pre = 4;
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

% data_inter = randintrlv(data_ldpc_codeRM,intrlv_state);
% data_inter = reshape(reshape(data_ldpc_codeRM, [], bits2sym).', 1,[]);
% sum(abs(data_inter-data_ldpc_codeRM))

% return
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
FAMP = 1100;
load mycustomfilter
tx = sdrtx('Pluto','RadioID','usb:0',filtnv{:});
tx.UseCustomFilter = true;
tx.filtCoefficients = [FAMP*32, zeros(1,31)];
tx.filtGain = 0;

tx.ShowAdvancedProperties = true;
tx.CenterFrequency    = fc; % Задаем несущую частоту
tx.BasebandSampleRate = BW; % Задаем полосу частот
tx.Gain = -15; %-89.75 to 0
transmitRepeat(tx, normalize(tx_chirp.'));    % Осуществляем непрерывную передачу


%% ================================= Прием с помощью Pluto
% figure(3)
for nIter=1:10

%         rx = sdrrx('Pluto','RadioID','usb:1');
        load mycustomfilter_rx
        rx = sdrrx('Pluto','RadioID','usb:1',filtnv_rx{:});
        rx.UseCustomFilter = true;
        rx.filtCoefficients = [FAMP*32, zeros(1,31)];
        rx.filtGain = 0;
% return
        % info(rx) %104473541196000508002e000e2d1ffaa3
        entered_offset = -2.5*round(BW/Base);
%         entered_offset = 0;
        rx.CenterFrequency    = fc + entered_offset; % Выбираем несущую частоту
        rx.BasebandSampleRate = BW; % Выбираем полосу частот
        rx.SamplesPerFrame = tx_length*10; % Принимаем в 10 раз больше чем отправили
        
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
        channel_chirp_sto = rx_data_in_time.';

        %% ================================= Correlation
        [cor,lags] = xcorr(channel_chirp_sto, downch);
        cor(round(end/2):end) = cor(round(end/2):end).*gausswin(length(cor(round(end/2):end))).';
        [max_amp, max_idx] = max(abs(cor));
        start = lags(max_idx);
        channel_chirp_corr = channel_chirp_sto(abs(start)+1:end);
%         channel_chirp_corr = channel_chirp_corr(Base+1:end);

        N = Base;
        for i=1:num_pre
            fourier_rx(i,:) = abs(fft( channel_chirp_corr(i*N-N+1+N:i*N+N).*downch ));
%             fourier(i*N-N+1:i*N) = abs(fft( channel_chirp_corr(i*N-N+1:i*N).*downch ));
        end

% figure(1)
% plot(real(tx_chirp))
% hold on
% plot(real(channel_chirp_corr))
%         return

        %% ================================= Frequency correction
%         channel_chirp_corr = channel_chirp_corr(Base+1:end);
%         rx_preamble = channel_chirp_corr(1:Base*num_pre);
%         corrected_signal = channel_chirp_corr(Base*num_pre+1:sig_length);
        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM(channel_chirp_corr, num_pre, sig_length);

%         figure(1)
%         
%         hold on
%         plot(real(mod_chirp))
%         plot(real(corrected_signal))
%         xlim([1 Base])
%         return

        N = Base;
        for i=1:num_pre
            fourier_fcor(i,:) = abs(fft( rx_preamble(i*N-N+1:i*N).*downch ));
        end
        for i=1:num_symRM
            fourier_data(i*N-N+1:i*N) = abs(fft( corrected_signal(i*N-N+1:i*N).*downch ));
        end

figure(1)
subplot(221)
stem(fourier_rx.')
title('rx')

subplot(222)
stem(fourier_fcor.')
title('fcor')

subplot(223)
stem(fourier_data.')
title('fdata')

subplot(224)
plot( (abs(fft(channel_chirp_corr))).')
title('spec')

        %% ================================= демодуляция
        aos = 3;
        [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( corrected_signal, num_symRM, tx_preamble, rx_preamble);
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
        est2 = freq_data{3};
        est3 = freq_data{4};
        est4 = freq_data{4};
        freq_shift_est = est1+est2+est3+est4;

        fprintf('Num Err: %d\n', bit_err)
        fprintf('STOint:  %d\n\n', STOint)

        fprintf('Offset:    %.2f\n', entered_offset)
        fprintf('Est 1:    %.2f\n', est1)
        fprintf('Est 2:    %.2f\n', est2)
        fprintf('Est 3:    %.2f\n', est3)
        fprintf('Est 4:    %.2f\n', est4)
        fprintf('Offset Est:  %.2f\n', freq_shift_est)
        fprintf('Offset Err:  %.2f\n\n', abs(freq_shift_est+entered_offset))

        if(bit_err>600)
            err_pos = find(check_data_no_gray~=sv_cor);
            fprintf('check G = %s\n', num2str(check_data(err_pos)))
            fprintf('sv G    = %s\n', num2str(sv(err_pos)))
            fprintf('check  = %s\n', num2str(check_data_no_gray(err_pos)))
            fprintf('sv_cor = %s\n', num2str(sv_cor(err_pos)))
            
            tx_chirp2 = tx_chirp(Base*(num_pre+1)+1:end);
            for i = 1:length(err_pos)
                fft_idl(i,:) = abs(fft(tx_chirp2(err_pos(i)*Base-Base+1: err_pos(i)*Base).*downch));
                fft_err(i,:) = abs(fft(corrected_signal(err_pos(i)*Base-Base+1: err_pos(i)*Base).*downch));
            
                sig_idl(i,:) = tx_chirp2(err_pos(i)*Base-Base+1: err_pos(i)*Base);
                sig_err(i,:) = corrected_signal(err_pos(i)*Base-Base+1: err_pos(i)*Base);
            end
            
            % c = 3;  % big foffset without correction
            c= 1;
% figure(1)
% subplot(221)
% stem(fourier_rx.')
% title('rx')
% 
% subplot(222)
% stem(fourier_fcor.')
% title('fcor')
% 
% subplot(223); hold on
% stem([zeros(1,Base*(err_pos(c)-1)), 0.5*(fft_err(c,:)/max(fft_err(c,:)))], 'LineWidth', 3)
% stem(fourier_data/max(fourier_data(c,:))); hold off
% title('fdata')
            
            figure(3)
            subplot(211)
            plot(real(normalize(sig_idl(c,:))))
            hold on
            plot(real(normalize(sig_err(c,:))))
            legend('ideal', 'err')
            
            subplot(212)
            stem( fft_idl(c,:)/max(fft_idl(c,:)), 'LineWidth', 3 )
            hold on
            stem( fft_err(c,:)/max(fft_err(c,:)) )
            % xlim([460 510])
            legend('ideal', 'err')

            release(tx);
            release(rx);
            return
        end
    % Estimate the BER for both methods
    BER(nIter) = bit_err;
end

release(tx);
release(rx);



save('rx_data_in_time.mat','channel_chirp_sto')
save('data.mat','data')

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