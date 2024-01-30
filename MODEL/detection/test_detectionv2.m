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
nIter = 100;
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
nIter = 1000;
Detect = 0;

tic

snr = [-16:0];

for n=1:8
    fprintf('Iter: %d\n', length(snr)-n+1) 

    tx_preamble = repmat(chirp, 1, n);
    tx_chirp = [tx_preamble, mod_chirp];
        

    Detect = 0;
    snr = -5;
    for iter=1:nIter
    
        % Delay
        delay = randi(Base*1);
        delay_full = length(mod_chirp)+delay;
        tx_chirp_shift = [zeros(1, delay), tx_chirp, zeros(1,Base*10)];
        
        % Noise
        noise_power = 20*log10(std(tx_chirp)./10^(snr/10));
        noise_vec = wgn(1, length(tx_chirp_shift), noise_power, 'complex');
        tx_chirp_noise = tx_chirp_shift+noise_vec;
        num_sym_tx = floor(length(tx_chirp_noise)./Base);
        
        % ~~~~~~~~ Preamble synchronization ~~~~~~~~
        delay_est = 0;
        fourier = 0;
        c = 1;
        for i = 1:num_sym_tx
%                 d = tx_chirp_noise(32*(i-1)+1:32*i+128-32).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
            d = tx_chirp_noise(Base*(i-1)+1:Base*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
%             fourier = abs(fft(d));            % переводим результат в область частот 
            fourier = fourier + abs(fft(d));
            
            if(c==n)
                [peak_code, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
                delay_est = Base-(indexMax-1);
                c = 1;
                if( delay_est==delay)
                    Detect = Detect+1;
                end
                break
%                 figure(1)
%                 plot(fourier)
%                 return
            end
        c=c+1;
%             peak_reg = [indexMax-1, peak_reg(1:end-1)];
%             decision = range(peak_reg);
%             if(decision<=1)
%                 if(peak_reg(1)==0)
%                     delay_est = i*Base-Base*7;
%                 else
%                     delay_est = i*Base-Base*7+(Base-peak_reg(1));
%                 end
%                 if( delay_est==delay)
%                     Detect = Detect+1;
%                 end
%                 break
%             end 
        end

    end

    % Estimate the BER for both methods
    BER(n) = Detect/nIter;
end
toc

% return
%%
n=1:8;
figure(1)
semilogy(n,BER,'-*','color','k');
% plot(n,10*log10(BER),'-*','color','k');
grid
ylabel('Вероятность обнаружения')
xlabel('Число передаваемых символов преамбулы')
hold off


% save('BER_PRE.mat','BER')

