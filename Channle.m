clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 2e6;
fc = 1000e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
Ts = LORA.Ts;

LORA2 = myLoRaClass(SF+1,BW);
sync_sym = LORA2.downch;

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
% data = zeros(1,numinfobits);
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
% sync_sym = resample(downch,10,1);
tx_chirp = [sync_sym, tx_downch, tx_preamble, mod_chirp, tx_preamble]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);




%% Channel

fps = BW/Base;
H = [1, 0.5, 0.6, 0.1];

freqshift = 5.4*fps;
dfi = freqshift*2*pi/BW;

for j = 1:tx_length
    fshifsin(j) = tx_chirp(j)*exp(1i*dfi*j);
end

Delay = randi([1 10]);
timeshift = [zeros(1,Delay),fshifsin,zeros(1,100)];

timeshift = conv(timeshift,H);


channel_chirp_sto = awgn(timeshift,0,'measured');
%         channel_chirp_sto = resample(channel_chirp_sto,1,resamp_factor);

        %% ================================= Correlation
        [channel_chirp_corr, cor] = LORA.CORRELATION(channel_chirp_sto, sync_sym, tx_length, 1);
        channel_chirp_corr = channel_chirp_corr(Base*2+1:end);


% figure(10)
% plot(abs(cor))
% return
        %% ================================= Frequency correction
        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM_v3(channel_chirp_corr, num_pre);
%         [corrected_signal, w] = MY_NLMS(Base/2, 0.1, tx_preamble, rx_preamble, corrected_signal);
        [corrected_signal] = LORA.STO_COMP(corrected_signal, num_pre);

        N = Base;
        rx_pre = channel_chirp_corr(pre_len+1:pre_len*2);
        for i=1:num_pre
            fourier_rx(i,:) = abs(fft( rx_pre(i*N-N+1:i*N).*downch ));
        end
        for i=1:num_pre
            fourier_fcor(i,:) = abs(fft( rx_preamble(i*N-N+1:i*N).*downch ));
        end
        for i=1:num_symRM
            fourier_data(i*N-N+1:i*N) = abs(fft( corrected_signal(i*N-N+1:i*N).*downch ));
        end
    
        figure(3)
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
        semilogy( (abs(fftshift(fft(channel_chirp_corr)))).')
        ylim([1 1000])
        title('spec')
        drawnow

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

%         fprintf('Num Err: %d\n', bit_err)
%         fprintf('STO:     %.2f\n\n', STO)
% 
%         fprintf('Offset:    %.2f\n', entered_offset)
%         fprintf('Est 1:    %.2f\n', est1)
%         fprintf('Est 2:    %.2f\n', est2)
%         fprintf('Offset Est:  %.2f\n', freq_shift_est)
%         fprintf('Offset Err:  %.2f\n\n', abs(freq_shift_est+entered_offset))

    % Estimate the BER for both methods
    BER = bit_err

% 
% release(tx);
% release(rx);

% profile viewer



