clc
clear all
close all

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 2e6;
snr = [-20:1:0];
nIter = 20;

% LORA = myLoRaClass_RSG(SF,BW); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LORA = myLoRaClass_true(SF,BW); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LORA.fir_win = 1;
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;
Ts = LORA.Ts;

num_pre = 4;
intrlv_state = 13;

% num_sym = 1000; 
% numinfobits = num_sym*bits2sym; 
% data = randi([0 1],1, numinfobits); 

%% ================================= LDPC coding
numinfobits = 1200;

% Message/Iformation bits
data = randi([0 1],1, numinfobits); 
num_sym = length(data)/rc;

%% ================================= Mодуляция
% RS
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified_crcrs(data, num_sym);

tx_preamble = repmat(chirp,1,num_pre);
tx_downch = repmat(downch,1,num_pre);
sync_sym = myLoRaClass_RSG(SF+1,BW).downch;

tx_chirp = [sync_sym, tx_downch, tx_preamble, mod_chirp];
tx_length = length(tx_chirp);


%% ================================= BER

% Channel
h11 = load('h11.mat').h;
% h = Channel(1,1,1);
% save('h11.mat','h')
% 
% figure(1)
% stem(abs(h))
% return
h11 = [h11 zeros(1, tx_length-length(h11))];
H11 = fft(h11);

tx_chirp_h = ifft( fft(tx_chirp).*H11 );
% tx_chirp_h = tx_chirp;

% вводим частотный сдвиг
fps = BW/Base;
freq_shift = fps*0.0;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

for j=1:tx_length
    tx_chirp_hf(j)=tx_chirp_h(j)*exp(1i*dphi*j);
end

% % Time delay
% delay = randi(50);
% tx_chirp_hft = [zeros(1,delay), tx_chirp_hf, zeros(1, tx_length)];
tx_chirp_hft = tx_chirp_hf;

tic
% snr = 10;

% c = 1;
% for i=-30:10
% [rxSig, var] = awgn(tx_chirp_hft, i, 'measured');
%     var_log(c) = 10*log10(var);
%     var_decision(c) = floor(var_log(c)*6/16);
%     c = c+1;
% end
% % if
% % var_decision(var_decision)
% % 
% g = gausswin(17, 1.5).';
% % g1 = g(6:12);
% g2 = g(3:15);
% g3 = g;
% stem( g )
% hold on
% stem( tukeywin(110, 10).')
% % stem( abs(fft(t(110, 11).')) )
% return

for n = 1:length(snr)
    fprintf('Iter left: %d\n', length(snr)-n+1) 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        rxSig = awgn(tx_chirp_hft, snr(n), 'measured');

        % Time sync
%         [rxSig_corr, cor] = LORA.CORRELATION(rxSig, sync_sym, tx_length);
%         rxSig_corr = rxSig_corr(Base*2+1:end);
        rxSig_corr = rxSig(Base*2+1:end);

        % Freq Sync
        rx_preamble = rxSig_corr(num_pre*Base+1:num_pre*2*Base);
        corrected_signal = rxSig_corr(num_pre*2*Base+1:end);
%         [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM_v3(rxSig_corr, num_pre);
%         corrected_signal = rxSig_corr(Base*num_pre*2+1:end);

        % Demodulation
        % RS
        [soft_bits, hard_bits, sv_rs, sv, fourier, fourier_rs] = LORA.delorax_crcrs( corrected_signal, num_sym, tx_preamble, rx_preamble);


        % подсчет БЕР с учетом задержки
        err = sum(hard_bits~=data);

        % Increment the error and bit counters
        numErr = numErr + err;        
        NumData = NumData + numinfobits;
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



% return
% 
% save('lora_rsg_mod_ber.mat','BER')
% save('lora_rsg_ber.mat','BER')
% save('lora_rs_ber.mat','BER')
% save('snr_crc.mat','snr')

