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

LORA = myLoRaClass_true(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;
ts = LORA.ts;

num_sym = 1;
numinfobits = SF*num_sym;
data = randi([0 1],1, numinfobits); 

%% ================================= Канал (AWGN + Phase shift)
fps = BW/Base;
max_peak_shift = 6;
resamp_factor = 10;

snr = 0;
num_pre_list = 2:2:8;

tic
for n = 1:length(num_pre_list)

    %% =================================  Modulation
%     n_pre = num_pre_list(n);
n_pre = 8;

    [mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data, num_sym, 1);
    tx_preamble = repmat(chirp,1,n_pre);
    tx_downch = repmat(downch,1,n_pre);
    
    tx_chirp = [tx_downch, tx_preamble, mod_chirp];
    tx_length = length(tx_chirp);

    tx_chirp_fshift = zeros(1,tx_length);
    est_err_list = zeros(1,nIter);

    for iter = 1:nIter
        freq_shift = randi([-round(fps*max_peak_shift), round(fps*max_peak_shift)]);
        freq_shift = fps*1;
        dphi = freq_shift*2*pi*ts;% сдвиг
        dphi = 0;
        
        % вводим частотный сдвиг
        for j=1:tx_length
            tx_chirp_fshift(j)=tx_chirp(j)*exp(1i*dphi*j);
        end


%         time_shift = round(-max_peak_shift + 2*max_peak_shift*rand(), 1);
        time_shift = 0.3
        tx_chirp_ftshift = resample(tx_chirp_fshift, resamp_factor, 1);
        tx_chirp_ftshift = circshift(tx_chirp_ftshift, time_shift*10);
        tx_chirp_ftshift_dec = resample(tx_chirp_ftshift, 1, resamp_factor);
%         tx_chirp_ftshift_dec = circshift(tx_chirp_fshift, time_shift);

%         [freq_data, ~, corrected_preamb] = LORA.LORA_FREQ_ESTIM_v3(tx_chirp_ftshift_dec, n_pre);

rx_preamb = tx_chirp_ftshift_dec(n_pre*Base+1:n_pre*Base*2);
rx_downch = tx_chirp_ftshift_dec(1:n_pre*Base);



% ~~~~~~~~ 0. Coarse CFO and STO estimation ~~~~~~~~ 
            N = Base;
%             fps = BW/N;
%             r1 = 0;
%             r2 = 0;
            OS = 10;
%             for i = 1:n_pre
%                 r1 = r1 + abs( [fft(rx_preamb(i*N-N+1:N*i).*downch), zeros(1,N*(OS-1))] );
%                 r2 = r2 + abs( [fft(rx_downch(i*N-N+1:N*i).* chirp), zeros(1,N*(OS-1))] );
%             end
% r1 = [fft(rx_preamb(1:Base).*downch), zeros(1,N*(OS-1))];%[fft(chirp.*downch), zeros(1,N*(OS-1))];
% r1 = [fft(chirp.*downch), zeros(1,N*(OS-1))];
r1 = [fft(rx_preamb(1:Base)), zeros(1,N*(OS-1))].*[fft(conj(rx_preamb(Base+1:2*Base))), zeros(1,N*(OS-1))];
[val, idx] = max(abs(ifft(r1)));
idx2 = ((idx-1))/20
% idx2 = (1280-(idx-1))/10
            
        figure(1); hold on
%         plot(real(tx_chirp))
%         plot(real(tx_chirp_ftshift_dec))
%         xlim([1 512])
        plot(abs(ifft(r1)))
%         plot(abs(r2))
        return

        
            [~, fup_idx] = max(abs(r1));
            [~, fdown_idx] = max(abs(r2));
            
            if(fup_idx>(Base*OS/2))
                fup_idx = (fup_idx-1)-Base*OS;
            else
                fup_idx = (fup_idx-1);
            end
            if(fdown_idx>(Base*OS/2))
                fdown_idx = (fdown_idx-1)-Base*OS;
            else
                fdown_idx = (fdown_idx-1);
            end
            fup_idx = fup_idx/OS; 
            fdown_idx = fdown_idx/OS;
            
            STO = (fup_idx-fdown_idx)/2;
            STOint = round(STO);
            STOfraq = STO-STOint;
            CFO = fps*(fup_idx+fdown_idx)/2;
            CFOdphi = CFO*2*pi*(1/BW);

            STO
            STOint-STOfraq
            
% return
% 
        figure(1); hold on
%         plot(real(tx_chirp_ftshift_dec(1:n_pre*Base)))
        plot(real(tx_chirp))
        plot(real(corrected_preamb))
        xlim([1 512])
        return
% 
%         figure(1); hold on
% %         plot(real(tx_chirp_ftshift_dec(1:n_pre*Base)))
%         plot(real(tx_chirp))
%         plot(real(corrected_preamb))
%         xlim([1 512])
%         return

        
        %% =================================  AWGN
        tx_chirp_ftshift_dec_n = awgn(tx_chirp_ftshift_dec, snr, 'measured');

        %% =================================  Freq estimation
        [freq_data, ~, ~] = LORA.LORA_FREQ_ESTIM_v3(tx_chirp_ftshift_dec_n, n_pre);
    
        STO_est = freq_data{1};
        est_full = freq_data{2}+freq_data{3};
        est1 = freq_data{2};
        est2 = freq_data{3};

%         fprintf('STO    :  %.2f\n', time_shift)
%         fprintf('STO Est:  %.2f\n', STO_est)
%         fprintf('STO Est Err:  %.2f\n\n', abs(STO_est+time_shift))
        est_err_list(iter) = freq_shift-est1-est2;
%         return

    end
    est_err_from_pre(n) = std(est_err_list);
end

toc   

figure(1)
histogram(est_err_list, 60)

% figure(2)
% stem(est_err_list)
% ylabel('Frequency, Hz');
% grid on

figure(3)
plot(num_pre_list, est_err_from_pre)
ylabel('Frequency, Hz');
grid on


% save('ftest_new_std.mat','est_err_from_pre'); % сохранение неразвернутой преамбулы



% %         time_shift
% rx_preamb = tx_chirp_ftshift_dec(n_pre*Base+1:n_pre*Base*2);
% rx_downch = tx_chirp_ftshift_dec(1:n_pre*Base);
% 
% 
% 
% % ~~~~~~~~ 0. Coarse CFO and STO estimation ~~~~~~~~ 
%             N = Base;
%             fps = BW/N;
%             r1 = 0;
%             r2 = 0;
%             OS = 8;
%             for i = 1:n_pre
%                 r1 = r1 + abs( fft([rx_preamb(i*N-N+1:N*i).*downch, zeros(1,N*(OS-1))]) );
%                 r2 = r2 + abs( fft([rx_downch(i*N-N+1:N*i).* chirp, zeros(1,N*(OS-1))]) );
%             end
%             
% %         figure(1); hold on
% %         plot(real(tx_chirp))
% %         plot(real(tx_chirp_ftshift_dec))
% %         xlim([1 512])
% % %         plot(abs(r1))
% % %         plot(abs(r2))
% %         return
%             [~, fup_idx] = max(abs(r1));
%             [~, fdown_idx] = max(abs(r2));
%             
%             if(fup_idx>(Base*OS/2))
%                 fup_idx = (fup_idx-1)-Base*OS;
%             else
%                 fup_idx = (fup_idx-1);
%             end
%             if(fdown_idx>(Base*OS/2))
%                 fdown_idx = (fdown_idx-1)-Base*OS;
%             else
%                 fdown_idx = (fdown_idx-1);
%             end
%             fup_idx = fup_idx/OS; 
%             fdown_idx = fdown_idx/OS;
%             
%             STO = (fup_idx-fdown_idx)/2;
%             STOint = round(STO);
%             STOfraq = STO-STOint;
%             CFO = fps*(fup_idx+fdown_idx)/2;
%             CFOdphi = CFO*2*pi*(1/BW);
% 
%             STO
%             STOint-STOfraq
%             
% return