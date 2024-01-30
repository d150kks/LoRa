clc
clear all
close all



%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 125e3;
nIter = 1000;
snr = 0;

LORA = myLoRaClass_true(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;
ts = LORA.ts;

num_sym = 1000;
numinfobits = SF*num_sym;
data = randi([0 1],1, numinfobits); 

%% =================================  Modulation
% n_pre = num_pre_list(n);
n_pre = 4;
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data, num_sym, 1);
% tx_preamble = repmat(chirp,1,n_pre);
% tx_downch = repmat(downch,1,n_pre);

LORA2 = myLoRaClass_true(SF+3,BW);
tx_preamble = LORA2.chirp;
% tx_preamble = chirp2
tx_chirp = [tx_preamble, mod_chirp];
signal_length = length(tx_chirp);
% tx_chirp = [tx_preamble, tx_preamble, circshift(chirp,64), mod_chirp];
% num_sym_tx = length(tx_chirp)./Base;



%% ================================= Detection
% num_sym_tx = floor(length(tx_chirp_noise)./Base);

% ~~~~~~~~ Demodulation ~~~~~~~~
nIter = 100;
Detect = 0;

tic

snr = [-16:0];
% snr=10;
for n=1:length(snr)
    fprintf('Iter: %d\n', length(snr)-n+1) 
    [numErr, NumData] = deal(0);
    
    for iter=1:nIter
    
        % Delay
        delay = randi(Base*1);
    %     delay=107;
        delay_full = length(mod_chirp)+delay;
        tx_chirp_shift = [zeros(1, delay), tx_chirp, zeros(1,length(tx_chirp))];
        
        % Noise
    %     snr = -0;
        noise_power = 20*log10(std(tx_chirp)./10^(snr(n)/10));
        noise_vec = wgn(1, length(tx_chirp_shift), noise_power, 'complex');
        tx_chirp_noise = tx_chirp_shift+noise_vec;
        num_sym_tx = floor(length(tx_chirp_noise)./Base);
%         ar = [];
        
        % ~~~~~~~~ Preamble synchronization ~~~~~~~~
        [output_signal, cor] = LORA.CORRELATION(tx_chirp_noise, tx_preamble, signal_length);
        rxSig = output_signal(length(tx_preamble)+1:end);

%         figure(1)
%         plot(real(mod_chirp))
%         hold on
%         plot(real(rxSig))
%         xlim([1 512])
%         return

    
        % ~~~~~~~~ Demodulation ~~~~~~~~
    %     rxSig = tx_chirp_noise(delay_est+1+Base*8:delay_est+Base*8+length(mod_chirp));
%         rxSig = detect_sig(sfd_start:sfd_start-1+length(mod_chirp));
        [soft_bits, hard_bits, sv_decode, sv, fourier] = LORA.delorax_modified( rxSig, num_sym);
    
        % подсчет БЕР с учетом задержки
        err = sum(hard_bits~=data);
    
        % Increment the error and bit counters
        numErr = numErr + err;        
        NumData = NumData + numinfobits;

%         if(numErr~=0)
%             figure(1)
%             plot(abs(ar))
% 
%             figure(2)
%             stem(abs(fft_cor))
%             break
%         end
    end

    % Estimate the BER for both methods
    BER(n) = numErr/NumData;
end
toc

% return
%%
figure(1)
semilogy(snr,BER,'-*','color','k');
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');


save('BER_SNRSYNCH.mat','BER')

return