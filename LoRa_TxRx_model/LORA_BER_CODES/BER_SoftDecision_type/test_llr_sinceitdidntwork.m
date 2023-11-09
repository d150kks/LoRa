clc
clear all
close all

% return
tic

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 30e6;
snr = [-16:1:0];
nIter = 20;

LORA = myLoRaClass_test(SF,BW);
Base = LORA.Base;
downch = LORA.downch;

num_pre = 8;
intrlv_state = 13;

% num_sym = 10;
% nbits = 1200; 
% data = randi([0 1],1, nbits); 
% data = [0 0 1];

grayCode = LORA.grayCode;
M0 = LORA.M0;
M1 = LORA.M1;
% return

%% ================================= LDPC coding
[cfgLDPCEnc,cfgLDPCDec] = LORA.generateConfigLDPC(3/4, 1944);
        
% Number of message bits
numinfobits = cfgLDPCEnc.NumInformationBits; 
numcodebits = cfgLDPCEnc.BlockLength; 

% Message/Iformation bits
data = randi([0 1],1, numinfobits); 
% data = zeros(1,numinfobits);
data_ldpc = ldpcEncode(data.', cfgLDPCEnc).';
% data_ldpc_intrlv = randintrlv(data_ldpc, intrlv_state);


%% ================================= Rate matching
[dataRM, numcodebitsRM, num_sym, zeros2end, flagRM] = LORA.RM(data_ldpc);


%% ================================= CRC coding
[data_code] = LORA.codeCRC(dataRM, num_sym);

%% ================================= Mодуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data_code, num_sym, 1);
tx_preamble = repmat(LORA.chirp,1,num_pre);

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
freq_shift = fps*0.0; %%%%%%%%%%%%%%%%%%%%%%%%
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(mod_chirp)
    mod_chirp(j)=mod_chirp(j)*exp(1i*dphi*j);
end


snr = 0;

% АБГШ 
rxSig = awgn(mod_chirp,snr,'measured');
[rx_preamble, nvar] = awgn(tx_preamble, snr, 'measured');

% демодуляция
hard_bits = zeros(1,SF*num_sym);
soft_bits = [];
sv = zeros(1,num_sym);
sv_cor = zeros(1,num_sym);

% ~~~~~~~~ Demodulation ~~~~~~~~
for i = 1:num_sym

    d = rxSig(Base*i-Base+1:Base*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
    fourier = abs(fft(d));            % переводим результат в область частот
    [peakMak, indexMax] = max( fourier );

    % ~~~~~~~~ Hard Decisions ~~~~~~~~
    [peak_notcor, peak_cor, peakMakcor, dbits] = LORA.HARD_CRC_DEMOD(fourier);
    sv(i) = peak_notcor;
    sv_cor(i) = peak_cor;
    hard_bits(SF*i-SF+1:SF*i) = dbits;


    % ~~~~~~~~ Soft Decisions ~~~~~~~~
    fourier_soft = fourier(grayCode+1);

    for nBit=1:SF
        m0 = M0(:, nBit);
        m0(m0==0)=[];
        m1 = M1(:, nBit);
        m1(m1==0)=[];
        LLR = -(1/nvar)*(min( (peakMak-fourier(grayCode(m0)+1)).^2 ) - min( (peakMak-fourier(grayCode(m1)+1)).^2 ));
%                 LLR = -(1/nvar)*(min( (peakMakcor-fourier_soft(m0)).^2 ) - min( (peakMakcor-fourier_soft(m1)).^2 ));
%                 LLR = -(1/nvar)*(min( (sv_cor(i)-m0).^2 ) - min( (sv_cor(i)-m1).^2 ));
%                 soft_bits(i*SF-SF+nBit) = LLR;
        soft_bits = [soft_bits, LLR];
    end

end
% peakMakcor
figure(2); hold on
% plot(fourier_soft)
plot( -normalize(soft_bits) )
plot( normalize(hard_bits) )
xlim([1 100])
return
        % декодирование CRC
        [rx_data_ldpc_intrlv] = LORA.decodeCRC(soft_bits, num_sym, zeros2end, flagRM);
        
        % LDPC decoding
        maxnumiter = 10;
        rx_data_ldpc = randdeintrlv(rx_data_ldpc_intrlv, intrlv_state);
        data_decode = ldpcDecode(rx_data_ldpc.', cfgLDPCDec, maxnumiter).';

        % подсчет БЕР с учетом задержки
        err = sum(data_decode~=data);


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
% save('lora_crc_softAMP_ber.mat','BER')
save('lora_crc_softPOS_ber.mat','BER')
% save('snr_crc.mat','snr')

