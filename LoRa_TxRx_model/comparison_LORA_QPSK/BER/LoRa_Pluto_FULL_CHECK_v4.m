clc
clear all
close all

tic

%% ================================= Переменные

% коэффициенты
SF = 11;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 30e6;
fc = 930e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
Ts = LORA.Ts;

intrlv_state = 13;
num_pre = 4;
warning('off')
% return

% profile on
%% ================================= LDPC coding
[cfgLDPCEnc,cfgLDPCDec] = LORA.generateConfigLDPC(3/4, 1944);
        
% Number of message bits
numinfobits = cfgLDPCEnc.NumInformationBits; 
numcodebits = cfgLDPCEnc.BlockLength; 

% Message/Iformation bits
data = randi([0 1],1, numinfobits); 
data = zeros(1,numinfobits);
data_ldpc = ldpcEncode(data.', cfgLDPCEnc).';
data_ldpc_intrlv = randintrlv(data_ldpc, intrlv_state);


%% ================================= Rate matching
[data_ldpc_intrlvRM, numcodebitsRM, num_symRM, zeros2end, flagRM] = LORA.RM(data_ldpc_intrlv);

%% ================================= CRC coding
[data_crc_ldpc_codeRM] = LORA.codeCRC(data_ldpc_intrlvRM, num_symRM);

%% ================================= Mодуляция
[mod_chirp, check_data, check_data_no_gray] = LORA.lorax_modified( data_crc_ldpc_codeRM, num_symRM, 1);
tx_preamble = repmat(LORA.chirp,1,num_pre);
pre_len = num_pre*Base;
tx_downch = repmat(downch,1,num_pre);
sync_sym = resample(downch,2,1);
tx_chirp = [sync_sym, tx_downch, tx_preamble, mod_chirp, tx_downch, tx_preamble]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);

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
tx.Gain = -20; %-89.75 to 0
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

        entered_offset = 2.5*round(BW/Base);
%         entered_offset = 0;
        rx.CenterFrequency    = fc + entered_offset; % Выбираем несущую частоту
        rx.BasebandSampleRate = BW+100; % Выбираем полосу частот
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
        [channel_chirp_corr, cor] = LORA.CORRELATION(channel_chirp_sto, sync_sym, tx_length);
        channel_chirp_corr = channel_chirp_corr(Base*2:end);

        N = Base;
        rx_pre = channel_chirp_corr(pre_len+1:pre_len*2);
        for i=1:num_pre
            fourier_rx(i,:) = abs(fft( rx_pre(i*N-N+1:i*N).*downch ));
        end

        %% ================================= Frequency correction
        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM_v3(channel_chirp_corr, num_pre);
        [corrected_signal] = LORA.STO_COMP(corrected_signal, num_pre);

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
        [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( corrected_signal, num_symRM, tx_preamble, rx_preamble);

        
        %% ================================= декодирование CRC
        [rx_data_ldpc_intrlv] = LORA.decodeCRC(soft_bits, num_symRM, zeros2end, flagRM);

        
        %% ================================= LDPC decoding
        maxnumiter = 100;
        rx_data_ldpc = randdeintrlv(rx_data_ldpc_intrlv, intrlv_state);
        rx_data = ldpcDecode(rx_data_ldpc.', cfgLDPCDec, maxnumiter).';


        %% ================================= Debugging
        bit_err = sum(abs(data~=rx_data));

        STO = freq_data{1};
        est1 = freq_data{2};
        est2 = freq_data{3};
        freq_shift_est = est1+est2;

        fprintf('Num Err: %d\n', bit_err)
        fprintf('STO:     %.2f\n\n', STO)

        fprintf('Offset:    %.2f\n', entered_offset)
        fprintf('Est 1:    %.2f\n', est1)
        fprintf('Est 2:    %.2f\n', est2)
        fprintf('Offset Est:  %.2f\n', freq_shift_est)
        fprintf('Offset Err:  %.2f\n\n', abs(freq_shift_est+entered_offset))

%                     release(tx);
%             release(rx);
%             return
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

% profile viewer

% save('data_from_pluto_null.mat','channel_chirp_sto')
save('data_from_pluto_100.mat','channel_chirp_sto')
save('data_pluto.mat','data')

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