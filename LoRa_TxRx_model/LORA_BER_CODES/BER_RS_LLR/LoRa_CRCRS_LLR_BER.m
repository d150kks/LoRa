clc
clear all
close all


%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
Base_rc = 2^rc;
rc_factor = 2^SF/Base_rc;
bits2sym = rc;
BW = 2;
snr = [-16:1:0];
nIter = 10;

LORA = myLoRaClass_true(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;

num_pre = 8;
intrlv_state = 13;

num_sym = 1000;
numinfobits = num_sym*bits2sym; 
data = randi([0 1],1, numinfobits); 

M0rs = LORA.M0rs;
M1rs = LORA.M1rs;
grayCode = LORA.grayCode;
RS = rc;
rs_factor = rc_factor;
rs_peaks_lut = (0:Base_rc-1);

% %% ================================= LDPC coding
% [cfgLDPCEnc,cfgLDPCDec] = LORA.generateConfigLDPC(3/4, 1944);
%         
% % Number of message bits
% numinfobits = cfgLDPCEnc.NumInformationBits; 
% numcodebits = cfgLDPCEnc.BlockLength; 
% 
% % Message/Iformation bits
% data = randi([0 1],1, numinfobits); 
% data_ldpc = ldpcEncode(data.', cfgLDPCEnc).';
% data_ldpc_intrlv = randintrlv(data_ldpc, intrlv_state);
% num_sym = length(data_ldpc_intrlv)/rc;

%% ================================= Mодуляция
% [mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data, num_sym, 1);
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified_crcrs(data, num_sym);
tx_preamble = repmat(chirp,1,num_pre);

%% ================================= BER
tic

h11 = zeros(1, length(mod_chirp));
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);

% mod_chirp = ifft( fft(mod_chirp).*H11 );
mod_chirp = mod_chirp;

fps = BW/Base;
freq_shift = fps*0.0;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(mod_chirp)
    mod_chirp(j)=mod_chirp(j)*exp(1i*dphi*j);
end

snr = 10;
for n = 1:length(snr)
%     fprintf('Iter: %d\n', n) 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        rxSig = awgn(mod_chirp,snr(n),'measured');
        [rx_preamble, var] = awgn(tx_preamble,snr(n),'measured');

        % демодуляция
        nvar = var;  
        
        hard_bits = zeros(1,RS*num_sym);
        soft_bits = zeros(1,RS*num_sym);
        sv = zeros(1,num_sym);
        sv_rs = zeros(1,num_sym);

        % Demodulation
        for i=1:num_sym
            
            d = rxSig(Base*i-Base+1:Base*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
            fourier = abs(fft(d));            % переводим результат в область частот
            fourier_rs = LORA.reduced_set_fourier(fourier);

            [peak_RS, indexMax] = max( fourier_rs ); % находим щелчок  частоты в чирпе
            sv(i) = rs_peaks_lut(grayCode==(indexMax-1));
%             sv(i) = grayCode(indexMax);
            sv_rs(i) = sv(i)*rs_factor;

            [~, indexMax2] = max( fourier_rs(grayCode+1) );
            fprintf('indexMax: %d\n', indexMax)
            fprintf('peredal: %d\n', check_no_gray(1))
            fprintf('sv(i): %d\n', sv(i)) 
            fprintf('indexMax2-1: %d\n', indexMax2-1) 
%             return

            % ~~~~~~~~ Soft and Hard Decisions ~~~~~~~~
            % Hard
            hard_bits(RS*i - RS+1:RS*i) = int2bit(sv(i).', RS).';

            % Soft Decisions
            for nBit=1:RS
                m0 = M0rs(:, nBit);
                m0(m0==0)=[];
                m0g = grayCode(m0)+1;

                m1 = M1rs(:, nBit);
                m1(m1==0)=[];
                m1g = grayCode(m1)+1;
                LLR = -(1/nvar)*(min( ( peak_RS-fourier_rs(m0g)).^2 ) - min( ( peak_RS-fourier_rs(m1g)).^2 ));
%                 LLR = -(1/nvar)*(min( ( peak_RS-fourier_rs(grayCode(m0)+1)).^2 ) - min( ( peak_RS-fourier_rs(grayCode(m1)+1)).^2 ));
%                 LLR = -(1/nvar)*(min( ( peak_RS-fourier_soft(m0)).^2 ) - min( ( peak_RS-fourier_soft(m1)).^2 ));
                soft_bits(i*RS-RS+nBit) = LLR;
            end
        end

figure(1)
subplot(211); hold on
% stem( gausswin(16+1).' )
stem( fourier_rs )
% stem(fourier_soft)

subplot(212); hold on
plot( -normalize(soft_bits) )
plot( normalize(hard_bits))
xlim([1 100])
return
%         [soft_bits, hard_bits, sv_rs, sv, fourier, fourier_rs] = LORA.delorax_crcrs( rxSig, num_sym);
%         [soft_bits, hard_bits, sv_rs, sv, fourier, fourier_rs] = LORA.delorax_crcrs( rxSig, num_sym, tx_preamble, rx_preamble);

%         % LDPC decoding
%         maxnumiter = 100;
%         rx_data_ldpc = randdeintrlv(soft_bits, intrlv_state);
%         data_decode = ldpcDecode(soft_bits.', cfgLDPCDec, maxnumiter).';



        % подсчет БЕР с учетом задержки
        err = sum(hard_bits~=data);
%         err = sum(data_decode~=data);

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
% save('lora_crcrs_ber.mat','BER')
% save('lora_rs_ber2.mat','BER')
% save('snr_crc.mat','snr')

