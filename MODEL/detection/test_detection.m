clc
clear all
close all


% upch*upch/upch*downch
% eight consequtive symbols
% large window

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 125e3;
nIter = 10000;
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
tx_preamble = repmat(chirp,1,n_pre);
tx_downch = repmat(downch,1,n_pre);

% SFD = [downch, downch, downch(1:32)];
SFD = [downch, downch];
% SFD = downch;
SFD_demapper = conj(SFD);

% tx_chirp = [tx_downch, tx_preamble, zeros(1,Base), mod_chirp];
tx_chirp = [tx_preamble, tx_preamble, SFD, mod_chirp];
% tx_chirp = [tx_preamble, tx_preamble, circshift(chirp,64), mod_chirp];
% num_sym_tx = length(tx_chirp)./Base;

%% ================================= Канал (AWGN + Phase shift)
% % Delay
% delay = randi(Base*1);
% % delay = 96;
% delay_full = length(mod_chirp)+delay;
% tx_chirp_shift = [mod_chirp, zeros(1, delay), tx_chirp];
% 
% Noise
% snr = -0;
% noise_power = 20*log10(std(tx_chirp)./10^(snr/10));
% noise_vec = wgn(1, length(tx_chirp_shift), noise_power, 'complex');
% % 10*log10(std(tx_chirp)/std(noise_vec))
% tx_chirp_noise = tx_chirp_shift+noise_vec;



%% ================================= Detection
% num_sym_tx = floor(length(tx_chirp_noise)./Base);

% ~~~~~~~~ Demodulation ~~~~~~~~
nIter = 100;
Detect = 0;

tic

snr = [-16:0];
% snr=0;
for n=1:length(snr)
    fprintf('Iter: %d\n', length(snr)-n+1) 
    [numErr, NumData] = deal(0);
    
    for iter=1:nIter
    
        % Delay
        delay = randi(Base*1);
    %     delay=107;
        delay_full = length(mod_chirp)+delay;
        tx_chirp_shift = [zeros(1, delay), tx_chirp, zeros(1,Base*10)];
        
        % Noise
    %     snr = -0;
        noise_power = 20*log10(std(tx_chirp)./10^(snr(n)/10));
        noise_vec = wgn(1, length(tx_chirp_shift), noise_power, 'complex');
        tx_chirp_noise = tx_chirp_shift+noise_vec;
        num_sym_tx = floor(length(tx_chirp_noise)./Base);
%         ar = [];
        
            % ~~~~~~~~ Preamble synchronization ~~~~~~~~
            peak_reg = [-0.0301    0.6508   -0.6244   -1.7215   -0.7358   -0.4550    0.9543]*100;
            delay_est = 0;
            fourier = 0;
%             c = 0;
            for i = 1:num_sym_tx
%                 d = tx_chirp_noise(32*(i-1)+1:32*i+128-32).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
                d = tx_chirp_noise(Base*(i-1)+1:Base*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
                fourier = abs(fft(d));            % переводим результат в область частот 
%                 fourier = fourier + abs(fft(d));

%                 ar = [ar fourier];
                [peak_code, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
            
                peak_reg = [indexMax-1, peak_reg(1:end-1)];
                decision = range(peak_reg);
                if(decision<=2)
                    if(peak_reg(1)==0)
                        delay_est = i*Base-Base*7;
                    else
                        delay_est = i*Base-Base*7+(Base-peak_reg(1));
                    end
                    if( delay_est==delay)
                        Detect = Detect+1;
                    end
                    break
                end 
            end
    
        
        % ~~~~~~~~ Start Frame Delimiter ~~~~~~~~
    %     detect_sig = tx_chirp_noise(delay_est+1+Base*7:delay_est+Base*7+length(mod_chirp)+length(SFD));
        detect_sig = tx_chirp_noise(delay_est+1+Base*7:end);
%         fft_cor = ifft(fft(detect_sig(1:Base*3)).*fft([chirp, zeros(1,Base*2)]));
%         fft_cor = ifft(fft(detect_sig(1:Base*5)).*fft([SFD_demapper, zeros(1,Base*3-32)]));
        fft_cor = ifft(fft(detect_sig(1:Base*4)).*fft([SFD_demapper, zeros(1,Base*2)]));
        [~, sfd_start] = max(abs(fft_cor));
    
%     figure(1)
%     stem(abs(fft_cor))
    % plot(imag(detect_sig))
    % 
    % return
    
        % ~~~~~~~~ Demodulation ~~~~~~~~
    %     rxSig = tx_chirp_noise(delay_est+1+Base*8:delay_est+Base*8+length(mod_chirp));
        rxSig = detect_sig(sfd_start:sfd_start-1+length(mod_chirp));
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


save('BER_FFTSYNCH.mat','BER')

return


% %% ================================= Detection
% % ~~~~~~~~ Demodulation ~~~~~~~~
% Detect = 0;
% nIter = 10000;
% preamble = [chirp, chirp, chirp, downch, downch, chirp, downch];
% % preamble = [chirp, chirp, chirp, -chirp, -chirp, chirp, -chirp];
% tx_chirp = [preamble, zeros(1,Base), mod_chirp];
% 
% tic
% for iter=1:nIter
% 
% % Delay
% delay = randi(Base*1);
% delay_full = length(mod_chirp)+delay;
% tx_chirp_shift = [mod_chirp, zeros(1, delay), tx_chirp];
% 
% % Noise
% snr = -6;
% noise_power = 20*log10(std(tx_chirp)./10^(snr/10));
% noise_vec = wgn(1, length(tx_chirp_shift), noise_power, 'complex');
% tx_chirp_noise = tx_chirp_shift+noise_vec;
% 
% 
% % ~~~~~~~~ Original synchronization ~~~~~~~~
% [cor,lags] = xcorr(tx_chirp_noise, preamble);
% [~, max_idx] = max(abs(cor));
% delay_est = lags(max_idx);
% if(delay_est==delay_full)
%     Detect = Detect+1;
% end
% end
% 
% TD = 100*Detect/nIter;
% toc
% 
% plot(abs(cor))
% % return



% %% ================================= Detection
% % ~~~~~~~~ Demodulation ~~~~~~~~
% Detect = 0;
% nIter = 10000;
% % preamble = [chirp, chirp, chirp, downch, downch, chirp, downch];
% % preamble = [chirp, chirp, chirp, -chirp, -chirp, chirp, -chirp];
% tx_chirp = [tx_preamble, tx_preamble, zeros(1,Base), mod_chirp];
% 
% % Delay
% delay = randi(Base*1);
% delay = 111;
% delay_full = length(mod_chirp)+delay;
% tx_chirp_shift = [mod_chirp, zeros(1, delay), tx_chirp];
% 
% % Noise
% snr = 10;
% noise_power = 20*log10(std(tx_chirp)./10^(snr/10));
% noise_vec = wgn(1, length(tx_chirp_shift), noise_power, 'complex');
% tx_chirp_noise = tx_chirp_shift+noise_vec;
% num_sym_tx = floor(length(tx_chirp_noise)./Base);
% 
% ar = [];
% peak_reg = zeros(1,8);
% for i = 1:num_sym_tx
%         d = tx_chirp_noise(Base*i-Base+1:Base*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
%         fourier = abs(fft(d));            % переводим результат в область частот %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ar = [ar fourier];
%         [peak_code, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
%     
%         % ~~~~~~~~ Original synchronization ~~~~~~~~
%         peak_reg = [indexMax-1, peak_reg(1:end-1)];
%         decision = range(peak_reg);
%         fprintf('peak register: %s\n', num2str(peak_reg))
%         if(decision==0)
%             if(peak_reg(1)==0)
%                 delay_est = i*Base-Base*8;
%             else
%                 delay_est = i*Base-Base*8+(Base-peak_reg(1));
%             end
%             if(delay_est==delay_full)
%                 Detect = Detect+1;
%             end
%         plot(abs(ar))
%         return
%         end 
% end
% 
% TD = 100*Detect/nIter;
% toc
% 
% plot(real(ar))
% return














% large_downch= repmat(downch, 1, 8);
% large_downch= [tx_preamble, tx_preamble];
% large_fft = xcorr(large_downch, tx_chirp);
% for i=1:110
%     large_demod = (large_downch.*tx_chirp_noise(1+i:8*Base+i));
%     large_fft(:,i) = fft(large_demod);
% end

xup = [0, sum(real(chirp))];
yup = [0, sum(imag(chirp))];
xd = [0, sum(real(downch))];
yd = [0, sum(imag(downch))];



h11 = zeros(1, Base);
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);

chirp_h = ifft(fft(chirp.*H11));

large_fft = fft(awgn(tx_preamble, -5, 'measured').*tx_downch);
mean(large_fft.*conj(large_fft))

 LORA2 = myLoRaClass_true(SF+1,BW);
 chirp2 =LORA2.chirp;

n=1:Base;
n2=1:Base*2;
% n=zeros(1,Base);
figure(3); hold on
% plot(abs(large_fft))
% plot(xup,yup,...
%      xd,yd)
plot3(n, real(chirp).', imag(chirp).', '-')
% plot3(n2, real(chirp2), imag(chirp2))
% plot(chirp_h)



return

points2check = [5:8];

figure(1)
subplot(211)
stem(sv)
text(points2check, sv(points2check), num2str(sv(points2check).') )
subplot(212)
plot(sv_amp)

% figure(2)
% plot(sv, sv_amp, '.')




% x = (1:10)' ;
% y = rand(size(x)) ;
% plot(x,y,'.r')
% hold on
% text(x,y,num2str(y))


% %% ================================= Detection my
% num_sym_tx = floor(length(tx_chirp_noise)./(Base))-8;
% large_downch = [tx_preamble, tx_downch];
% % large_downch = repmat(downch, 1, 8);
% 
% % ~~~~~~~~ Demodulation ~~~~~~~~
% counter = 1;
% % fourier_array = zeros(1,num_sym_tx*8*Base);
% fourier_array = [];
% for i = 1:num_sym_tx
% 
%     d = tx_chirp_noise(Base*(i-1)+1: Base*(8+i-1)).*large_downch;   % перемножаем входной и опорный ОБРАТНый чирп
%     fourier = abs(fft(d));            % переводим результат в область частот %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fourier_array = [fourier_array, fourier];
% %     fourier_array(8*Base*i-8*Base+1:8*Base*i) = fourier;
% %     [peak_code, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
% 
%     % ~~~~~~~~ Synchronization ~~~~~~~~
% 
% 
%     % ~~~~~~~~ Debugging ~~~~~~~~
% %     clc
% %     fprintf('decision: %d\n', decision)
% %     fprintf('peak register: %s\n', num2str(peak_reg))
% %     if(mod(i,8)==0)
% end
% 
% plot(abs(fourier_array))
% [~, delay_est1] = max(fourier_array)
% return