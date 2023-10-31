classdef myLoRaClass
   
    properties
        % коэффициенты
        SF;         % коэффициент расширения спектра (от 7 до 12)
        Base;       % База сигнала
        BW;         % Signal Bandwidth
        Ts;         % длительность сигнала

        % время и частота
        ts;         % время дискретизации
        m;          % коэффициент измнения частоты
        t;          % полоса времени ОДНОГО чирпа

        % Создание чирпов
        chirp;      % upchirp
        downch;     % downchirp

        % Gray code
        grayCode;   % Gray Code

        % CRC
        crclen = 4; % length of crc code
        bits2sym;   % bits per symbol
        aos = 3;    % Area of Peak search to Calculate CRC
        aos_win;    % Window of Peak search to Calculate CRC

        %
        M0;         % matrix with peak-to-bitvalue==0
        M1;         % matrix with peak-to-bitvalue==1
    end

    methods
        %% Constructor
        function obj = myLoRaClass(SF, BW)

            % ~~~~~~~~ Chirp parameters ~~~~~~~~
            obj.SF = SF; 
            obj.Base = 2^SF;
            obj.BW = BW; 
            obj.Ts = (2^SF)/BW;
            obj.ts = 1/(BW);              % время дискретизации
            obj.m = BW/obj.Ts;              % коэффициент измнения частоты
            obj.t = 0:obj.ts:obj.Ts-obj.ts;         % полоса времени ОДНОГО чирпа

            % ~~~~~~~~ Ref Chirps ~~~~~~~~
            obj.chirp  = exp( 1i * (2*pi*(0 + (obj.m*obj.t.^2)/2)) );
            obj.downch = exp(-1i * (2*pi*(0 + (obj.m*obj.t.^2)/2)) );

            % ~~~~~~~~ Gray code ~~~~~~~~
            obj.grayCode = zeros(1, obj.Base);
            for i = 0:obj.Base-1
                obj.grayCode(i+1) = bitxor(i, bitshift(i, -1));
            end

            % ~~~~~~~~ CRC code ~~~~~~~~
            obj.bits2sym = obj.SF-obj.crclen;
            obj.aos_win  = -obj.aos:obj.aos;

            % ~~~~~~~~ LLR Map ~~~~~~~~
            bitmap = de2bi(0:obj.Base-1);
            obj.M0 = zeros(obj.Base, obj.SF);
            obj.M1 = zeros(obj.Base, obj.SF);
            for nBit=1:obj.SF
                for nSym = 1:obj.Base
                    if( bitmap(nSym, nBit) == 0)
                        obj.M0(nSym,nBit)= nSym;
                    else
                        obj.M1(nSym,nBit)= nSym;
                    end
                end
            end
            
        end


        %% ================================= Rate matching
        function [dataRM, numbitsRM, num_symRM, zeros2end, flagRM] = RM(obj, data)
            numbits = length(data);
            num_sym = numbits/obj.bits2sym;
            bits2end = mod( numbits, obj.bits2sym);
            zeros2end = obj.bits2sym-bits2end;

            if( bits2end~=0 )
                numbitsRM = numbits+zeros2end;
                dataRM = [data, zeros(1,zeros2end)];
                num_symRM = numbitsRM/obj.bits2sym;
                flagRM = 1;
            else
                flagRM = 0;
                numbitsRM = numbits;
                dataRM = data;
                num_symRM = num_sym;
            end
        end


        %% ================================= CRC coding
        % ------------------------- CRC Coding -------------------------
        function [data_code] = codeCRC(obj, data, num_sym)

            % CRC coding
            data_code = zeros(1, num_sym*obj.SF);
            for i=1:num_sym
                data_win = data(i*obj.bits2sym-obj.bits2sym+1:i*obj.bits2sym);
                CRC_Bits = obj.CRC4(data_win.').';
                data_code(i*obj.SF-obj.SF+1:i*obj.SF) = [data_win, CRC_Bits].';
            end

        end

        %% ================================= CRC decoding
        % ------------------------- CRC Decoding -------------------------
        function [data_decode] = decodeCRC(obj, data, num_sym, zeros2end, flagRM)

            data_decode = zeros(1, num_sym*obj.bits2sym);
            for i=1:num_sym
                data_win = data(i*obj.SF-obj.SF+1:i*obj.SF);
                data_decode(i*obj.bits2sym-obj.bits2sym+1:i*obj.bits2sym) = data_win(1:obj.bits2sym);
            end
            if(flagRM==1)
                data_decode = data_decode(1:end-zeros2end);
            end

        end

        %% ================================= lorax_modified
        function [mod_chirp, check_data, check_no_gray] = lorax_modified(obj, data, num_sym, gray)

            % формирование пустых массивов
            check_data = zeros(1,num_sym);
            check_no_gray = zeros(1,num_sym);
            mod_chirp = zeros(1, num_sym*obj.Base);           % чирпы
            
            % Модуляция
            for i = 1:num_sym
                code_word = bi2de(data(obj.SF*i-obj.SF+1:obj.SF*i));      % значение кодового слова
                if(gray==1)
                    code_word_gray = obj.grayCode(code_word+1);
                else
                    code_word_gray = code_word;
                end
                check_data(i)    = code_word_gray; % значение закодированного символа
                check_no_gray(i) = code_word;
                cs = single((code_word_gray/obj.Base)*obj.Ts/obj.ts);             % место сдвига
                
                chirp1 = obj.chirp(cs+1:end); % первый чирп
                chirp2 = obj.chirp(1:cs); % второй чирп

                mod_chirp(i*obj.Base-obj.Base+1:obj.Base*i) = [chirp1, chirp2];
            end 
        end

        %% ================================= delorax_modified
        function [sv_decode, sv, fourier] = delorax_modified(obj, mod_chirp, num_sym)

            % Demodulation
            sv = zeros(1,num_sym);
            sv_decode = zeros(1,num_sym);
            for i = 1:num_sym
                d = mod_chirp(obj.Base*i-obj.Base+1:obj.Base*i).*obj.downch;   % перемножаем входной и опорный ОБРАТНый чирп
                
                fourier = abs(fft(d));            % переводим результат в область частот
                [~, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
                [~, indexMaxGray] = max( fourier(obj.grayCode+1) ); % находим щелчок  частоты в чирпе
    
                % вычисляем значение кодового слова исходя из базы сигнала
                sv(i) = indexMax-1;
                sv_decode(i) = indexMaxGray-1;
            end
            
        end
   

        %% ================================= DELORAX_CRC
        function [sv, sv_cor, peakMakcor, dbits] = HARD_CRC_DEMOD(obj, fourier)

            % ~~~~~~~~ Initial conditions ~~~~~~~~
            check_crc = 1;

            % ~~~~~~~~ CRC Demodulation ~~~~~~~~
            while check_crc~=0
            
                % Circular shift in peak search
                [~, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
                sv = indexMax;
                peak_win = indexMax+obj.aos_win;
                for pk=1:length(peak_win)
                    if(peak_win(pk)<=0)
                        peak_win(pk)=peak_win(pk)+obj.Base;
                    end
                    if(peak_win(pk)>obj.Base)
                        peak_win(pk)=peak_win(pk)-obj.Base;
                    end
                end

                % Set type of the peak sort
                peaks_amp = fourier(peak_win);
                peak_sort = zeros(1,2*obj.aos+1);
                peak_sort(1) = indexMax;
                for pk=1:obj.aos
                    if( peaks_amp(obj.aos+1+pk)>peaks_amp(obj.aos+1-pk) )
                        peak_sort(2*pk)       = peak_win(obj.aos+1+pk);
                        peak_sort(2*(pk+1)-1) = peak_win(obj.aos+1-pk);
                    else
                        peak_sort(2*pk)       = peak_win(obj.aos+1-pk);
                        peak_sort(2*(pk+1)-1) = peak_win(obj.aos+1+pk);
                    end
                end
                sort_amp = fourier(peak_sort);
                peak_sort(sort_amp==0)=[];
                sort_amp(sort_amp==0)=[];
    
                % Search peak according to the crc
                for n=1:length(peak_sort)
                    sv_cor = find(obj.grayCode==peak_sort(n)-1)-1;
                    dbits = de2bi(sv_cor,obj.SF)';
                    dbits = dbits(:)';
                    check_crc = sum(obj.CRC4(dbits.').');

                    if(check_crc==0)
                        peakMakcor = sort_amp(n);
                        break
                    else
                        fourier(peak_sort(n))=0;
                    end
                end
            end
        end

        % demodulation 
        function [soft_bits, hard_bits, sv, sv_cor, fourier] = DELORAX_CRC(obj, mod_chirp, num_sym, tx_preamble, rx_preamble)

            % ~~~~~~~~ Params ~~~~~~~~
            nvar = std( abs(tx_preamble-rx_preamble).^2);
            hard_bits = zeros(1,obj.SF*num_sym);
            soft_bits = zeros(1,obj.SF*num_sym);
            sv = zeros(1,num_sym);
            sv_cor = zeros(1,num_sym);

            % ~~~~~~~~ Demodulation ~~~~~~~~
            for i = 1:num_sym
            
                d = mod_chirp(obj.Base*i-obj.Base+1:obj.Base*i).*obj.downch;   % перемножаем входной и опорный ОБРАТНый чирп
                fourier = abs(fft(d));            % переводим результат в область частот

                % ~~~~~~~~ Hard Decisions ~~~~~~~~
                [peak_notcor, peak_cor, peakMakcor, dbits] = obj.HARD_CRC_DEMOD(fourier);
                sv(i) = peak_notcor;
                sv_cor(i) = peak_cor;
                hard_bits(obj.SF*i-obj.SF+1:obj.SF*i) = dbits;
            
                % ~~~~~~~~ Soft Decisions ~~~~~~~~
                fourier_soft = fourier(obj.grayCode+1);
                for nBit=1:obj.SF
                    m0 = obj.M0(:, nBit);
                    m0(m0==0)=[];
                    m1 = obj.M1(:, nBit);
                    m1(m1==0)=[];
                    LLR = -(1/nvar)*(min( (peakMakcor-fourier_soft(m0)).^2 ) - min( (peakMakcor-fourier_soft(m1)).^2 ));
                    soft_bits(i*obj.SF-obj.SF+nBit) = LLR;
                end

            end
        end

        %% ================================= LORA_FREQ_ESTIM
        function [r1, r2, fup_idx, fdown_idx, phi1, phi2] = fraq(obj, rx_preamb, rx_downch)
            % 
            N = obj.Base;
            OS = 10;
            r1 = 0;
            num_pre = 4;
            for i = 1:num_pre
                r1 = r1 + abs( fft([rx_preamb(i*N-N+1:N*i).*obj.downch, zeros(1,N*(OS-1))]) );
            end

            r2 = fft( [rx_downch.*conj(obj.downch), zeros(1,obj.Base*(OS-1))] );
            [~, fup_idx] = max(abs(r1));
            [~, fdown_idx] = max(abs(r2));

            if(fup_idx>(obj.Base*OS/2))
                fup_idx = (fup_idx-1)-obj.Base*OS;
            end
            if(fdown_idx>(obj.Base*OS/2))
                fdown_idx = (fdown_idx-1)-obj.Base*OS;
            end
%             fup_idx = round(fup_idx/OS); 
%             fdown_idx = round(fdown_idx/OS);
            fup_idx = fup_idx/OS; 
            fdown_idx = fdown_idx/OS; 
            %###################

            phi2=0;
            phi1=0;
        end

        function [STOint, STOfraq, CFO, CFOdphi] = fraq2(obj, rx_preamb, rx_downch, num_pre)
            % ~~~~~~~~ 0. Coarse CFO and STO estimation ~~~~~~~~ 
            N = obj.Base;
            fps = obj.BW/N;
            OS = 10;
            r1 = 0;
            r2 = 0;
            for i = 1:num_pre
                r1 = r1 + abs( fft([rx_preamb(i*N-N+1:N*i).*obj.downch, zeros(1,N*(OS-1))]) );
                r2 = r2 + abs( fft([rx_downch(i*N-N+1:N*i).* obj.chirp, zeros(1,N*(OS-1))]) );
            end
            
%             r2 = fft( [rx_downch.*conj(obj.downch), zeros(1,obj.Base*(OS-1))] );
            [~, fup_idx] = max(abs(r1));
            [~, fdown_idx] = max(abs(r2));
            
            if(fup_idx>(obj.Base*OS/2))
                fup_idx = (fup_idx-1)-obj.Base*OS;
            end
            if(fdown_idx>(obj.Base*OS/2))
                fdown_idx = (fdown_idx-1)-obj.Base*OS;
            end
            fup_idx = fup_idx/OS; 
            fdown_idx = fdown_idx/OS;
            
            STO = (fup_idx-fdown_idx)/2;
            STOint = round(STO);
            STOfraq = STO-STOint;
            CFO = fps*(fup_idx+fdown_idx)/2;
            CFOdphi = CFO*2*pi*(1/obj.BW);
        end

        %% ================================= Correlation
        function [output_signal, cor] = CORRELATION(obj, input_signal, preamble, signal_length, conv_or_fft)

            switch conv_or_fft
                case 1
                    [cor,lags] = xcorr(input_signal, preamble);
                    cor(round(end/2):end) = cor(round(end/2):end).*gausswin(length(cor(round(end/2):end)), 1).';
                    [~, max_idx] = max(abs(cor));
                    figure
                    plot(abs(cor))
                    start = lags(max_idx);
                    output_signal = input_signal(start+1:start+signal_length);
                    
                case 0 
                    if( log2(signal_length)<ceil(log2(signal_length)) || log2(signal_length)<10)
                        padding = 2^ceil(log2(signal_length))-signal_length;
                    else
                        padding=0;
                    end
                    input_signal_fft = fft([input_signal, zeros(1,padding)]);
                    preamble_fft = fft(preamble, length(input_signal_fft));
                    cor = ifft(input_signal_fft.*conj(preamble_fft));
                    cor = cor.*gausswin(length(input_signal_fft)).';
                    length(cor)
                    [~, start] = max(abs(cor));

                    figure
                    plot(abs(cor))

                    output_signal = input_signal(start:start+signal_length-1);

                otherwise
                    fprintf("ASS")
            end
        end



        %% ================================= Dphi Compensation
        function [output_signal] = DPHI_COMP(obj, input_signal, dphi)
            len = length(input_signal);
            output_signal = zeros(1, len);
            for j=1:len
                output_signal(j)=input_signal(j).*exp(1i*dphi*j*(-1));
            end
        end

        %% ================================= SAMPLING ERROR EST AND COMP
        function [output_signal] = STO_COMP(obj, input_signal, num_pre)

            % input params
            N = obj.Base;
            fps = obj.BW/obj.Base;
            num_symRM = length(input_signal)/obj.Base;
            preamble_end = input_signal(end-num_pre*N*2+1:end);
            downch_end = preamble_end(1:num_pre*N);
            preamb_end = preamble_end(num_pre*N+1:num_pre*N*2);
% num_symRM


            [STOint, STOfraq, CFO, CFOdphi] = obj.fraq2(preamb_end, downch_end, num_pre);
            STO = STOint+STOfraq;
            est_te = STO*fps+CFO;
            
%             [est_te, ~] = obj.COARSE_FREQ_ESTIM(preamble_end, num_pre);
            dphi_te = est_te*2*pi*(1/obj.BW)/num_symRM; % сдвиг
            dphi_pdt = dphi_te;
            
            output_signal = zeros(1,length(input_signal));
            for j=1:length(input_signal)
                output_signal(j) = input_signal(j).*exp(-1i*j*(dphi_pdt));
                if(mod(j,N)==0)
                    dphi_pdt = dphi_pdt+dphi_te;
                end
            end
            output_signal = output_signal(1:end-num_pre*N);
        end

        %% ================================= CYC SHIFT
        function [output_vector] = CYC_SHIFT(obj, input_vector)
            len = length(input_vector);
            output_vector = zeros(1,len);
            for j=1:len
                if(input_vector(j)>obj.Base)
                    output_vector(j)=input_vector(j)-obj.Base;
                elseif(input_vector(j)<=0)
                    output_vector(j)=input_vector(j)+obj.Base;
                else
                    output_vector(j)=input_vector(j);
                end
            end
        end

        %% ================================= Coarse freq estimation
        function [est1, dphi1] = COARSE_FREQ_ESTIM(obj, input_signal, num_pre)
            % ~~~~~~~~ 1. Coarse estimation ~~~~~~~~ 

            N = obj.Base;
            fps = obj.BW/N;
            OS = 5;
            NOS = N*OS;
            fourier = 0;

            for i = 1:num_pre
                fourier = fourier + fftshift( abs(fft( [input_signal(i*N-N+1:N*i).*obj.downch, zeros(1,N*(OS-1))] )) );
            end
            [~, ind1] = max( fourier );
            pre_align = (ind1-1)-NOS/2;  

            est1 = pre_align*fps/OS; % UNCOMMENT
            dphi1 = est1*2*pi*(1/obj.BW); % сдвиг
        end

        %% ================================= 4-x Stage Freq Est and Comp

        function [freq_data, corrected_signal, corrected_preamb] = LORA_FREQ_ESTIM_v3(obj, input_signal, num_pre)
        
            % ~~~~~~~~ Description ~~~~~~~~ 
            %   input_signal - preamble+payload signals
            %   N - length of the one chirp
            %   num_pre - num of the preambles
            
            % ~~~~~~~~ Initializtion ~~~~~~~~
            % Extract payload and preamble signals
            N = obj.Base;
            fps = obj.BW/obj.Base;
            pre_len = obj.Base*num_pre;
            rx_control = input_signal(1:pre_len*2);
            rx_downch = rx_control(1:pre_len);
            rx_preamb = rx_control(pre_len+1:pre_len*2);

            % ~~~~~~~~ 0. Coarse CFO and STO estimation ~~~~~~~~ 
            [STOint, STOfraq, CFO, CFOdphi] = obj.fraq2( rx_preamb, rx_downch, num_pre);
%             STOint, STOfraq, CFO, CFOdphi
            STO = STOint+STOfraq;

            % ~~~~~~~~ 00. Coarse CFO and STO compensation ~~~~~~~~ 
            % ~~~~~~~~ 01. STO compensation ~~~~~~~~ 
            input_signal_fft = fft(input_signal);
            STOintphi = STO*2*pi/length(input_signal_fft);
%             STOintphi = STOint*2*pi/length(input_signal_fft);
            input_signal_ifft = ifft(obj.DPHI_COMP(input_signal_fft, STOintphi));

            % ~~~~~~~~ 02. STO compensation ~~~~~~~~ 
%             rx_preamb2 = input_signal_ifft(1:N*(num_pre+1));
%             rx_preamb2 = obj.DPHI_COMP(rx_preamb2, CFOdphi);

            rx_control2 = input_signal_ifft(1:pre_len*2);
            rx_control2 = obj.DPHI_COMP(rx_control2, CFOdphi);
            rx_preamb2 = rx_control2(pre_len+1:pre_len*2);

            
            
            
            
            % ~~~~~~~~ 3.0 Fine estimation ~~~~~~~~ 
            % Устранение фазового сдвига
            Yaos = -2:2;
            argumon = zeros(1, num_pre-1);
            for i = 1:num_pre-1
                Y1 = fft(rx_preamb2(i*N+1:N*i+N).*obj.downch);
                Y0 = fft(rx_preamb2(i*N-N+1:N*i).*obj.downch);
                [~, indmax1] = max(Y1);
                [~, indmax0] = max(Y0);
                indwin1 = obj.CYC_SHIFT(indmax1+Yaos);
                indwin0 = obj.CYC_SHIFT(indmax0+Yaos);
            
                argumon(i) = angle( sum( Y1(indwin1) .* conj(Y0(indwin0)) ) );
            end
            
            [arg] = mean(argumon);
            est3 = arg/(2*pi*obj.Ts);
            dphi3 = est3*2*pi/obj.BW; % сдвиг
            
%             % ~~~~~~~~ 3.1 Fine Compensation ~~~~~~~~ 
%             rx_preamb3 = obj.DPHI_COMP(rx_preamb2, dphi3);
% 
%             % ~~~~~~~~ 4. Coarse estimation X2 ~~~~~~~~ 
%             [est4, dphi4] = obj.COARSE_FREQ_ESTIM(rx_preamb3, num_pre);


            % ~~~~~~~~ 5.0 Final Compensation ~~~~~~~~ 
%             CFOdphi = 0;
%             dphi3 = 0;
            CFOdphi+dphi3;
            input_signal2 = obj.DPHI_COMP(input_signal_ifft, (CFOdphi+dphi3+0));
            corrected_preamb = input_signal2(pre_len+1:pre_len*2);
            corrected_signal = input_signal2(pre_len*2+1:end);
%             input_signal2 = obj.DPHI_COMP(input_signal_ifft, (CFOdphi+dphi3+0));
%             corrected_preamb = input_signal2(obj.Base+1:obj.Base*(num_pre+1));
%             corrected_signal = input_signal2(obj.Base*(num_pre+1)+1:end);

            % ~~~~~~~~ Debugging ~~~~~~~~ 
            freq_data = {STO, CFO, est3, 0};
        
        end

        %% ================================= CRC
        function [CRC_Bits] = CRC4(obj, input);
            % Тип CRC
            crc_poly = [1, 0, 0, 1, 1];    
            crc_len  = 4;
            % elseif(crc_type==3)
            %     crc_poly = [1, 0, 1, 1];    
            %     crc_len  = 3;
            % end
            crc_rem   = zeros(1,crc_len+1);
            tmp_array = vertcat(input, zeros(crc_len,1));
            
            % Вычисление CRC
            for n=1:length(input)+crc_len
                for k=1:crc_len 
                    crc_rem(k) = crc_rem(k+1); 
                end
                crc_rem(crc_len+1) = tmp_array(n);
            
                if(crc_rem(1) ~= 0)
                    for k=1:crc_len+1
                        crc_rem(k) = mod(crc_rem(k)+crc_poly(k), 2); 
                    end
                end
            end
            CRC_Bits = zeros(crc_len,1);
            for n=1:crc_len 
                CRC_Bits(n,1) = crc_rem(n+1); 
            end
        end

        function [TEraw] = TIME_ERR(obj, product);
            
            [~, ind_max] = max( abs(product) );
            if((ind_max-1)<=0)
                mag1 = abs(product(obj.Base+(ind_max-1)));
            else
                mag1 = abs(product(ind_max-1));
            end
        
            mag2 = abs(product(ind_max));
        %     mag3 = abs(product(ind_max+1));
            if((ind_max+1)>obj.Base)
                mag3 = abs(product((ind_max+1)-obj.Base));
            else
                mag3 = abs(product(ind_max+1));
            end
            TEraw = (mag3-mag1)/mag2;
        end

        %% ================================= LDPC
        function [cfgLDPCEnc,decodercfg] = generateConfigLDPC(obj, rate,varargin)
            
            % Checking input arguments and default inputs
            defaultCodewordLen = 648;
            defaultStandard = 'wlan';
            defaultDecoderAlgo = 'bp';
            
            expectedStandards = {'wlan'};
            expectedDecoderAlgos = {'bp','layered-bp','norm-min-sum','offset-min-sum'};
            
            parsVar = inputParser;
            
            errorMsgRate = 'Value must be positive and less than 1'; 
            validateRate= @(x) assert((x > 0) && (x <= 1),errorMsgRate);
            
            addRequired(parsVar,'rate',validateRate);
            
            errorMsgCWL = 'Value must be positive and less than 1'; 
            validateCW= @(x) assert((x > 0),errorMsgCWL);
            addOptional(parsVar,'codewordLen',defaultCodewordLen,validateCW);
            
            addParameter(parsVar,'standard',defaultStandard,@(x) any(validatestring(x,expectedStandards)));
            addParameter(parsVar,'decoderAlgo',defaultDecoderAlgo,@(x) any(validatestring(x,expectedDecoderAlgos)));
            
            parse(parsVar,rate,varargin{:});
            

            % ------------------------- "WLAN"--------------------------------
            if strcmpi(parsVar.Results.standard,'wlan')==1
            
                % Checking supported metrics
            
                % Check supported code-rate
                if rate~= 1/2 && rate~= 2/3 && rate~= 3/4 && rate~= 5/6
                    error('Unsupported Code-Rate');
                end
            
                % Check supported codeword block length
                wlanCodewordLens = [648,1296,1944];
            
                if ismember(parsVar.Results.codewordLen,wlanCodewordLens)==0
                        error('Unsupported Codeword Block Length for WLAN standard');
                else
                    n = parsVar.Results.codewordLen;
                end
            
                % Checking the codeword length condition.
                switch n
            
                    case 648
                        blockSize  = 27; % Variabel "Z" is used to denote in the standard
            
                         % Checking the codeword-rate condition.
                        switch rate
                            case 1/2
                                P = [ 0 -1 -1 -1  0  0 -1 -1  0 -1 -1  0  1  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
                                    22  0 -1 -1 17 -1  0  0 12 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1 -1
                                    6 -1  0 -1 10 -1 -1 -1 24 -1  0 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1
                                    2 -1 -1  0 20 -1 -1 -1 25  0 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1
                                    23 -1 -1 -1  3 -1 -1 -1  0 -1  9 11 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1
                                    24 -1 23  1 17 -1  3 -1 10 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1
                                    25 -1 -1 -1  8 -1 -1 -1  7 18 -1 -1  0 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1
                                    13 24 -1 -1  0 -1  8 -1  6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1
                                    7 20 -1 16 22 10 -1 -1 23 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1
                                    11 -1 -1 -1 19 -1 -1 -1 13 -1  3 17 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1
                                    25 -1  8 -1 23 18 -1 14  9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0
                                    3 -1 -1 -1 16 -1 -1  2 25  5 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0
                                    ];
            
            
                            case 2/3
                                P = [ 25 26 14 -1 20 -1  2 -1  4 -1 -1  8 -1 16 -1 18  1  0 -1 -1 -1 -1 -1 -1
                                    10  9 15 11 -1  0 -1  1 -1 -1 18 -1  8 -1 10 -1 -1  0  0 -1 -1 -1 -1 -1
                                    16  2 20 26 21 -1  6 -1  1 26 -1  7 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1
                                    10 13  5  0 -1  3 -1  7 -1 -1 26 -1 -1 13 -1 16 -1 -1 -1  0  0 -1 -1 -1
                                    23 14 24 -1 12 -1 19 -1 17 -1 -1 -1 20 -1 21 -1  0 -1 -1 -1  0  0 -1 -1
                                    6 22  9 20 -1 25 -1 17 -1  8 -1 14 -1 18 -1 -1 -1 -1 -1 -1 -1  0  0 -1
                                    14 23 21 11 20 -1 24 -1 18 -1 19 -1 -1 -1 -1 22 -1 -1 -1 -1 -1 -1  0  0
                                    17 11 11 20 -1 21 -1 26 -1  3 -1 -1 18 -1 26 -1  1 -1 -1 -1 -1 -1 -1  0
                                    ];
            
                            case 3/4
                                P = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
                                    25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
                                    25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
                                    9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
                                    24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
                                    2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
                                    ];
            
                            case 5/6
                                P = [17 13  8 21  9  3 18 12 10  0  4 15 19  2  5 10 26 19 13 13  1  0 -1 -1
                                    3 12 11 14 11 25  5 18  0  9  2 26 26 10 24  7 14 20  4  2 -1  0  0 -1
                                    22 16  4  3 10 21 12  5 21 14 19  5 -1  8  5 18 11  5  5 15  0 -1  0  0
                                    7  7 14 14  4 16 16 24 24 10  1  7 15  6 10 26  8 18 21 14  1 -1 -1  0
                                    ];
            
                        end
            
                    case 1296
                        blockSize  = 54; % Variabel "Z" is used to denote in the standard
            
                        switch rate
                            case 1/2
                                P = [40 -1 -1 -1 22 -1 49 23 43 -1 -1 -1  1  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
                                    50  1 -1 -1 48 35 -1 -1 13 -1 30 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1 -1
                                    39 50 -1 -1  4 -1  2 -1 -1 -1 -1 49 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1
                                    33 -1 -1 38 37 -1 -1  4  1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1
                                    45 -1 -1 -1  0 22 -1 -1 20 42 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1
                                    51 -1 -1 48 35 -1 -1 -1 44 -1 18 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1
                                    47 11 -1 -1 -1 17 -1 -1 51 -1 -1 -1  0 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1
                                    5 -1 25 -1  6 -1 45 -1 13 40 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1
                                    33 -1 -1 34 24 -1 -1 -1 23 -1 -1 46 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1
                                    1 -1 27 -1  1 -1 -1 -1 38 -1 44 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1
                                    -1 18 -1 -1 23 -1 -1  8  0 35 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0
                                    49 -1 17 -1 30 -1 -1 -1 34 -1 -1 19  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0
                                    ];
            
                            case 2/3
                                P = [39 31 22 43 -1 40  4 -1 11 -1 -1 50 -1 -1 -1  6  1  0 -1 -1 -1 -1 -1 -1
                                    25 52 41  2  6 -1 14 -1 34 -1 -1 -1 24 -1 37 -1 -1  0  0 -1 -1 -1 -1 -1
                                    43 31 29  0 21 -1 28 -1 -1  2 -1 -1  7 -1 17 -1 -1 -1  0  0 -1 -1 -1 -1
                                    20 33 48 -1  4 13 -1 26 -1 -1 22 -1 -1 46 42 -1 -1 -1 -1  0  0 -1 -1 -1
                                    45  7 18 51 12 25 -1 -1 -1 50 -1 -1  5 -1 -1 -1  0 -1 -1 -1  0  0 -1 -1
                                    35 40 32 16  5 -1 -1 18 -1 -1 43 51 -1 32 -1 -1 -1 -1 -1 -1 -1  0  0 -1
                                    9 24 13 22 28 -1 -1 37 -1 -1 25 -1 -1 52 -1 13 -1 -1 -1 -1 -1 -1  0  0
                                    32 22  4 21 16 -1 -1 -1 27 28 -1 38 -1 -1 -1  8  1 -1 -1 -1 -1 -1 -1  0
                                    ];
            
                            case 3/4
                                P = [39 40 51 41  3 29  8 36 -1 14 -1  6 -1 33 -1 11 -1  4  1  0 -1 -1 -1 -1
                                    48 21 47  9 48 35 51 -1 38 -1 28 -1 34 -1 50 -1 50 -1 -1  0  0 -1 -1 -1
                                    30 39 28 42 50 39  5 17 -1  6 -1 18 -1 20 -1 15 -1 40 -1 -1  0  0 -1 -1
                                    29  0  1 43 36 30 47 -1 49 -1 47 -1  3 -1 35 -1 34 -1  0 -1 -1  0  0 -1
                                    1 32 11 23 10 44 12  7 -1 48 -1  4 -1  9 -1 17 -1 16 -1 -1 -1 -1  0  0
                                    13  7 15 47 23 16 47 -1 43 -1 29 -1 52 -1  2 -1 53 -1  1 -1 -1 -1 -1  0
                                    ];
            
                            case 5/6
                                P = [48 29 37 52  2 16  6 14 53 31 34  5 18 42 53 31 45 -1 46 52  1  0 -1 -1
                                    17  4 30  7 43 11 24  6 14 21  6 39 17 40 47  7 15 41 19 -1 -1  0  0 -1
                                    7  2 51 31 46 23 16 11 53 40 10  7 46 53 33 35 -1 25 35 38  0 -1  0  0
                                    19 48 41  1 10  7 36 47  5 29 52 52 31 10 26  6  3  2 -1 51  1 -1 -1  0
                                    ];
            
                        end
            
                    case 1944
                        blockSize  = 81; % Variabel "Z" is used to denote in the standard
            
                        switch rate
                            case 1/2
                                P = [57 -1 -1 -1 50 -1 11 -1 50 -1 79 -1  1  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
                                    3 -1 28 -1  0 -1 -1 -1 55  7 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1 -1
                                    30 -1 -1 -1 24 37 -1 -1 56 14 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1
                                    62 53 -1 -1 53 -1 -1  3 35 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1
                                    40 -1 -1 20 66 -1 -1 22 28 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1
                                    0 -1 -1 -1  8 -1 42 -1 50 -1 -1  8 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1
                                    69 79 79 -1 -1 -1 56 -1 52 -1 -1 -1  0 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1
                                    65 -1 -1 -1 38 57 -1 -1 72 -1 27 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1
                                    64 -1 -1 -1 14 52 -1 -1 30 -1 -1 32 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1
                                    -1 45 -1 70  0 -1 -1 -1 77  9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1
                                    2 56 -1 57 35 -1 -1 -1 -1 -1 12 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0
                                    24 -1 61 -1 60 -1 -1 27 51 -1 -1 16  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0
                                    ];
            
                            case 2/3
                                P = [61 75  4 63 56 -1 -1 -1 -1 -1 -1  8 -1  2 17 25  1  0 -1 -1 -1 -1 -1 -1
                                    56 74 77 20 -1 -1 -1 64 24  4 67 -1  7 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1
                                    28 21 68 10  7 14 65 -1 -1 -1 23 -1 -1 -1 75 -1 -1 -1  0  0 -1 -1 -1 -1
                                    48 38 43 78 76 -1 -1 -1 -1  5 36 -1 15 72 -1 -1 -1 -1 -1  0  0 -1 -1 -1
                                    40  2 53 25 -1 52 62 -1 20 -1 -1 44 -1 -1 -1 -1  0 -1 -1 -1  0  0 -1 -1
                                    69 23 64 10 22 -1 21 -1 -1 -1 -1 -1 68 23 29 -1 -1 -1 -1 -1 -1  0  0 -1
                                    12  0 68 20 55 61 -1 40 -1 -1 -1 52 -1 -1 -1 44 -1 -1 -1 -1 -1 -1  0  0
                                    58  8 34 64 78 -1 -1 11 78 24 -1 -1 -1 -1 -1 58  1 -1 -1 -1 -1 -1 -1  0
                                    ];
            
                            case 3/4
                                P = [48 29 28 39  9 61 -1 -1 -1 63 45 80 -1 -1 -1 37 32 22  1  0 -1 -1 -1 -1
                                    4 49 42 48 11 30 -1 -1 -1 49 17 41 37 15 -1 54 -1 -1 -1  0  0 -1 -1 -1
                                    35 76 78 51 37 35 21 -1 17 64 -1 -1 -1 59  7 -1 -1 32 -1 -1  0  0 -1 -1
                                    9 65 44  9 54 56 73 34 42 -1 -1 -1 35 -1 -1 -1 46 39  0 -1 -1  0  0 -1
                                    3 62  7 80 68 26 -1 80 55 -1 36 -1 26 -1  9 -1 72 -1 -1 -1 -1 -1  0  0
                                    26 75 33 21 69 59  3 38 -1 -1 -1 35 -1 62 36 26 -1 -1  1 -1 -1 -1 -1  0
                                    ];
            
                            case 5/6
                                P = [13 48 80 66  4 74  7 30 76 52 37 60 -1 49 73 31 74 73 23 -1  1  0 -1 -1
                                    69 63 74 56 64 77 57 65  6 16 51 -1 64 -1 68  9 48 62 54 27 -1  0  0 -1
                                    51 15  0 80 24 25 42 54 44 71 71  9 67 35 -1 58 -1 29 -1 53  0 -1  0  0
                                    16 29 36 41 44 56 59 37 50 24 -1 65  4 65 52 -1  4 -1 73 52  1 -1 -1  0
                                    ];
                        end
            
                end
            
                pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);
                cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);
                decodercfg = ldpcDecoderConfig(pcmatrix,parsVar.Results.decoderAlgo);
            
            end
        end



    end
end