%%LORA_FREQ_ESTIM
function [freq_data, corrected_signal] = LORA_FREQ_ESTIM(input_signal, downch, N, num_pre, BW, Base, Ts)

    %% Description
    %   input_signal - preamble+payload signals
    %   N - length of the one chirp
    %   num_pre - num of the preambles
    
    %% Main body
    %% Extract payload and preamble signals
    channel_chirp_treshold = input_signal(1:num_pre*N);
    payload_signal = input_signal(num_pre*N+1:end);
    
    %% ================================= 1. Coarse estimation
    % Определяем Hz/samp
    fps = BW/Base;
    ts=1/BW;
    
    for i = 1:num_pre-1
        [max_a1, ind1] = max(fft(channel_chirp_treshold(i*N-N+1:N*i).*downch));
        pre_align(i) = ind1;
        
        if (pre_align(i)>N/2)
            pre_align(i) = pre_align(i)-N;
        end
    end
    
    est1 = (mean(pre_align-1))*fps;
    dphi1 = est1*2*pi*ts; % сдвиг
    
    for j=1:length(channel_chirp_treshold)
        channel_chirp_realign(j)=channel_chirp_treshold(j).*exp(1i*dphi1*j*(-1));
    end
    
    %% ================================= 2. Fraq estimation
    % По двум последовательным одинаковым чирпам как в (1.) вычисляем CFO_fraq
    left_half = channel_chirp_realign(1:N/2);
    left_ref = downch(1:N/2);
    
    right_half = channel_chirp_realign(N/2+1:N);
    right_ref = downch(N/2+1:N);
    
    bpf3 = fft(left_half.*left_ref);
    bpf4 = fft(right_half.*right_ref);
    [max_a3] = max(bpf3);
    [max_a4] = max(bpf4);
    
    a11=max(max_a3);
    a12=max(max_a4);
    
    est2 = (angle(a12)-angle(a11))/(pi*Ts);
    dphi2 = est2*2*pi*Ts/Base;
    
    % точное устранение фазового набега
    for j=1:length(channel_chirp_realign)
        channel_chirp_frac_est(j) = channel_chirp_realign(j).*exp(1i*dphi2*j*(-1));
    end
    
    %% ================================= 3. fine estimation
    % Устранение фазового сдвига
    for i = 1:num_pre-1
        argumon(i*N-N+1:N*i) = channel_chirp_frac_est(i*N-N+1:N*i).*conj(channel_chirp_frac_est(i*N+1:N*i+N));
    end
    
    [arg] = max(sum(argumon));
    est3 = -angle(arg)/(2*pi*Ts);
    
    % Устраняем CFO_fraq сигнала
    dphi3 = est3*2*pi*Ts/Base; % сдвиг
    
    
    %% ================================= Debugging
    est_full = est1+est2+est3;
    %%  freq_err = freq_shift-est_full;
    freq_data = {est1, est2, est3, est_full};
    
    %% ================================= Correcting Payload signal
    dphi_full = est_full*2*pi*Ts/Base; % сдвиг
    for j = 1:length(payload_signal)
        idx = num_pre*N+j;
        corrected_signal(j) = payload_signal(j).*exp(1i*dphi_full*idx*(-1));
    end

end
