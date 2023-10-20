clc
clear all
close all

return
tic

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 30e6;
snr = [-16:1:0];
nIter = 10;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;

num_pre = 8;

% num_sym = 10;
% nbits = 1200; 
% data = randi([0 1],1, nbits); 
data = [0 0 1];

%% ================================= Rate matching
[dataRM, numcodebitsRM, num_sym, zeros2end, flagRM] = LORA.RM(data);


%% ================================= CRC coding
[data_code] = LORA.codeCRC(dataRM, num_sym);

%% ================================= Mодуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data_code, num_sym, 1);
tx_preamble = repmat(LORA.chirp,1,num_pre);
tx_sig = [tx_preamble, mod_chirp];

%% ================================= BER

tic

h11 = zeros(1, length(mod_chirp));
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);

mod_chirp = ifft( fft(mod_chirp).*H11 );

fps = BW/Base;
freq_shift = fps*2.5; %%%%%%%%%%%%%%%%%%%%%%%%
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(mod_chirp)
    mod_chirp(j)=mod_chirp(j)*exp(1i*dphi*j);
end


snr = 10;
for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        rxSig = awgn(mod_chirp,snr(n),'measured');

        % демодуляция
        aos = 3;

        % Param
        Gcode = LORA.grayCode;
        hard_bits = zeros(1,SF*num_sym);

        aos_win = -aos:aos;
        sv = zeros(1,num_sym);
        sv_cor = zeros(1,num_sym);
        soft_bits = [];
        
        % Demodulation
        for i = 1:num_sym
        
            d = rxSig(Base*i-Base+1:Base*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп     
            fourier = abs(fft(d));            % переводим результат в область частот

            % Initial conditions
            check_crc = 1;
            while check_crc~=0
        
                [~, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
                sv(i) = indexMax;
                peak_win = indexMax+aos_win;
            
                % cyclic shift
                for pk=1:length(peak_win)
                    if(peak_win(pk)<=0)
                        peak_win(pk)=peak_win(pk)+Base;
                    end
                    if(peak_win(pk)>Base)
                        peak_win(pk)=peak_win(pk)-Base;
                    end
                end
            
                % Sort amps and peaks
                peaks_amp = fourier(peak_win);
                peak_sort = zeros(1,2*aos+1);
                peak_sort(1) = indexMax;
                for pk=1:aos
                    if( peaks_amp(aos+1+pk)>peaks_amp(aos+1-pk) )
                        peak_sort(2*pk)       = peak_win(aos+1+pk);
                        peak_sort(2*(pk+1)-1) = peak_win(aos+1-pk);
%                         peak_sort = [peak_sort, peak_win(aos+1+pk), peak_win(aos+1-pk)];
                    else
                        peak_sort(2*pk)       = peak_win(aos+1-pk);
                        peak_sort(2*(pk+1)-1) = peak_win(aos+1+pk);
%                         peak_sort = [peak_sort, peak_win(aos+1-pk), peak_win(aos+1+pk)];
                    end
                end
                sort_amp = fourier(peak_sort);

                for k=1:length(peak_sort)
                    sv_cor(i) = find(Gcode==peak_sort(k)-1)-1;
                    dbits = de2bi(sv_cor(i),SF)';
                    dbits = dbits(:)';
                    check_crc = sum(LORA.CRC4(dbits.').');
    %                 peakMakcor = sort_amp(k);


        fourier = abs(fft(rxSig.*downch));
        figure(1)
        plot(fourier)
        return
                    if(check_crc==0)
                        break
                    else
                        fourier(peak_sort(k))=0;
                    end
                end

            end
        hard_bits(SF*i-SF+1:SF*i) = dbits;
%         hard_bits = [hard_bits, dbits];
        end
        [data_decode] = LORA.decodeCRC(hard_bits, num_sym, zeros2end, flagRM);


        % подсчет БЕР с учетом задержки
        err = sum(data_decode~=data);

        % Increment the error and bit counters
        numErr = numErr + err;        
        NumData = NumData + nbits;
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
% save('lora_crc_ber.mat','BER')
% save('snr_crc.mat','snr')

