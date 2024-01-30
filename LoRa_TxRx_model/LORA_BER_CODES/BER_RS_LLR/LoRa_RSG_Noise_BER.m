clc
clear all
close all

return
%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 2e6;
snr = [-16:1:0];
nIter = 100;

LORA = myLoRaClass_true(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;

num_pre = 4;

num_sym = 300;
% numinfobits = num_sym*rc; 
numinfobits = num_sym*SF; 
data = randi([0 1],1, numinfobits); 


%% ================================= Mодуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified_crcrs(data, num_sym);
% [mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data, num_sym, 1);
tx_preamble = repmat(chirp, 1, num_pre);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tx_preamble = repmat(chirp,1,num_pre);
tx_downch = repmat(downch,1,num_pre);
sync_sym = myLoRaClass_true(SF+1,BW).downch;
tx_chirp = [sync_sym, tx_downch, tx_preamble, mod_chirp];

% save('tx_bits.mat','data')
% save('tx_chirp.mat','tx_chirp')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:num_sym
    mod_win(Base*(i-1)+1:Base*i) = mod_chirp(Base*(i-1)+1:Base*i).*triang(Base).';
end


figure
% plot(abs(fft(mod_win)))
plot(real(chirp))

output = papr(tx_chirp)

function output = papr(input)
    peake = max(abs(real(input)).^2);
    energy = sum(abs(real(input)).^2)/length(input);
    output = 10*log10(peake./energy);
end

% return
% 
% 
% %% ================================= BER
% % return
% tic
% % snr = 5;
% % snr = -3;
% for n = 1:length(snr)
%     fprintf('Iter: %d\n', n) 
% 
% %     sig_power = std(mod_chirp).^2;
% %     noise_power = sig_power/(10^(snr(n)/10));
% %     noise = wgn(1, length(mod_chirp), 10*log10(noise_power), 'complex');
% %     snr_est = 10*log10( std(mod_chirp).^2/std(noise).^2 )
% %     10*log10(noise_power)
% 
%     [numErr, NumData] = deal(0);
% % return
%     for iter = 1:nIter
%         
%         % АБГШ 
%         rxSig = awgn(mod_chirp, snr(n), 'measured');
% %         noise = wgn(1, length(mod_chirp), 10*log10(noise_power), 'complex');
% %         rxSig = mod_chirp + noise;
%         rx_preamble = awgn(tx_preamble,snr(n),'measured');
% 
%         % демодуляция
%         [soft_bits, hard_bits, sv_decode, sv, fourier] = LORA.delorax_modified( rxSig, num_sym, tx_preamble, rx_preamble);
% %         [soft_bits, hard_bits, sv_rs, sv, fourier, fourier_rs] = LORA.delorax_crcrs( rxSig, num_sym, tx_preamble, rx_preamble);
% 
% 
% %         figure(1)
% %         plot( normalize(hard_bits))
% %         hold on
% %         plot( -normalize(soft_bits))
% %         xlim([1 100])
% %         return
% 
%         % подсчет БЕР с учетом задержки
%         err = sum(hard_bits~=data);
% 
%         % Increment the error and bit counters
%         numErr = numErr + err;        
%         NumData = NumData + numinfobits;
%     end
% 
%     % Estimate the BER for both methods
%     BER(n) = numErr/NumData;
% end
% toc
% 
% %%
% figure(1)
% semilogy(snr,BER,'-*','color','k');
% grid
% xlabel('SNR (dB)')
% ylabel('Bit Error Rate')
% hold off
% title('SNR');
% 
% 
% % save('lora_ber.mat','BER')
% % save('lora_rsg_ber.mat','BER')