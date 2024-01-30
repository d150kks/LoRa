clc
clear all
close all


%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 2e6;
snr = [-6:1:10];
nIter = 100;


% rs_array = [1, 2, 3, 4];
for rs = 1:3:4
    LORA = myLoRaClass_true(SF,BW);
    LORA.rs_size = rs;
    LORA.fir_win = 1;
    Base = LORA.Base;
    downch = LORA.downch;
    chirp = LORA.chirp;
    
    num_pre = 4;
    
    num_sym = 1000;
    numinfobits = num_sym*rc; 
    data = randi([0 1],1, numinfobits); 
    
    
    %% ================================= Mодуляция
    [mod_chirp, check_data, check_no_gray] = LORA.lorax_modified_crcrs(data, num_sym);
    % [mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data, num_sym, 1);
    tx_preamble = repmat(chirp, 1, num_pre);

    h11 = zeros(1, length(mod_chirp));
    h11(1) = 1;
    h11(2) = 0.9;
    h11(3) = 0.8;
    h11(5) = 0.8;
    h11(7) = 0.7;
    h11(9) = 0.5;
    h11(11) = 0.45;
    h11(13) = 0.3;
    h11(15) = 0.2;
    h11(16) = 0.1;
    H11 = fft(h11);
    
    tx_chirp_h = ifft( fft(mod_chirp).*H11 );
    
    % вводим частотный сдвиг
    fps = BW/Base;
    freq_shift = fps*1.5;
    dphi=freq_shift*2*pi*(1/BW);% сдвиг
    
    for j=1:length(tx_chirp_h)
        tx_chirp_h(j)=tx_chirp_h(j)*exp(1i*dphi*j);
    end
    
    %% ================================= BER
    tic
    % snr = 5;
    for n = 1:length(snr)
        fprintf('Iter: %d\n', n) 
    
        [numErr, NumData] = deal(0);
    
        for iter = 1:nIter
            
            % АБГШ 
            rxSig = awgn(tx_chirp_h, snr(n), 'measured');
            rx_preamble = awgn(tx_preamble,snr(n),'measured');
    
            % демодуляция
    %         [soft_bits, hard_bits, sv_decode, sv, fourier] = LORA.delorax_modified( rxSig, num_sym, tx_preamble, rx_preamble);
            [soft_bits, hard_bits, sv_rs, sv, fourier, fourier_rs] = LORA.delorax_crcrs( rxSig, num_sym, tx_preamble, rx_preamble);
    
    
            % подсчет БЕР с учетом задержки
            err = sum(hard_bits~=data);
    
            % Increment the error and bit counters
            numErr = numErr + err;        
            NumData = NumData + numinfobits;
        end
    
        % Estimate the BER for both methods
        BER(rs, n) = numErr/NumData;
    end
    toc
end

%%
% figure(1)
% semilogy(snr,BER,'-*','color','k');
% grid
% xlabel('SNR (dB)')
% ylabel('Bit Error Rate')
% hold off
% title('SNR');

figure(1)
semilogy(snr,BER);
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');


% save('lora_ber.mat','BER')
% save('lora_rsg_ber.mat','BER')