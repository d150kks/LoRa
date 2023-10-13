function [Out, param] = Channel(In,BW,SNR,TimeShift,FreqShift)
    
    if (nargin < 2)
        BW = 125e3;
    end
    if (nargin < 3)
        SNR = 0; 
    end
    if (nargin < 4) 
        TimeShift = 10;
    end
    if (nargin < 5) 
        FreqShift = 1324;
    end
    %% Входные параметры
        param.SNR       = SNR;                                                          % Отношение сигнал/шум
        param.TimeShift = TimeShift;                                                    % Временной сдвиг
        param.FreqShift = FreqShift;                                                         % Частотный сдвиг
        param.NumOfRay  = 10;                                                           % Количество лучей в канале
        param.Ray       = [1, 0.9, 0.8, 0.8, 0.7, 0.5, 0.45, 0.3, 0.2, 0.1]*0.5;            % Затухание лучей
        param.Delay     = [0,   4,   6,   7,   8,  10,   12,  15,  16,  20]+TimeShift;  % Задержка лучей

    %% Частотный сдвиг

        dfi = param.FreqShift*2*pi/BW;
        In_FreqShift = zeros(length(In),1);
    
        for j = 1:length(In)
            In_FreqShift(j) = In(j)*exp(1i*dfi*j);
        end

    %% Временной сдвиг + Многолучевость

        In_TimeShift = zeros(length(In_FreqShift)+max(param.Delay),1);
        In_TimeShift(max(param.Delay)+1:end) = In_FreqShift;
    
        if(param.NumOfRay > 1)
            for i=2:param.Delay
                In_TimeShift = In_TimeShift+([zeros(param.Delay(i),1);In_FreqShift;zeros(length(In_TimeShift)-param.Delay(i)-length(In_FreqShift),1)]).*param.Ray(i);
            end
        end
    
        Out = In_TimeShift;
end