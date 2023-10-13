function GEN(Frame_time, VISA, fd)
CH1_source = 'OUTPUT1 ON'; % ВКЛ ВЫКЛ 1 канал
CH2_source = 'OUTPUT2 ON'; % ВКЛ ВЫКЛ 2 канал
CH1_filter = 'OFF'; %% ВКЛ ВЫКЛ фильтр на 1 канале
CH2_filter = 'OFF'; %% ВКЛ ВЫКЛ фильтр на 2 канале
offset1 = 0; % смещение 1 канала
offset2 = 0; % смещение 2 канала
load1 = 50; % ВХ сопротивление 1 канала
load2 = 50; % Вх сопротивление 2 канала
amp = 0.1; % амплитуда в вольтах, НЕ рекомендуется больше 2!!!

%This function connects to a 33522 and uploads an I and Q baseband signal
%opens and creates a visa session for communication with function generator
fgen = instrfind('Type', 'visa-tcpip', 'RsrcName', VISA, 'Tag', '');
if isempty(fgen)
    fgen = visa('AGILENT', VISA);
%     fgen = visadev(VISA);
else
    fclose(fgen);
    fgen = fgen(1);
end
fgen.Timeout = 1000;
get(fgen,{'OutputBufferSize'})

% %calculate output buffer size
buffer = length(Frame_time)*8;
set (fgen,'OutputBufferSize',(buffer+125));

get(fgen,{'OutputBufferSize'})

%open connection to 33522
try
   fopen(fgen);
catch exception %problem occurred throw error message
    uiwait(msgbox('Error occurred trying to connect to the 33522, verify correct IP address','Error Message','error'));
    rethrow(exception);
end

%Query Idendity string and report
fprintf (fgen, '*IDN?');
idn = fscanf (fgen);
fprintf (idn)
fprintf ('\n\n')
mes = ['Connected to ' idn ' sending waveforms.....'];
h = waitbar(0,mes);
%Reset instrument
fprintf (fgen, '*RST');
%turn it to column vector
I = single(real(Frame_time).');
Q = single(imag(Frame_time).');

% signal_in_time = single(signal_in_time');

%Clear volatile memory
fprintf(fgen, 'SOURce1:DATA:VOLatile:CLEar');
fprintf(fgen, 'SOURce2:DATA:VOLatile:CLEar');

% fprintf(fgen, 'CALibration:SECure:STATe:ON'); %% !!SYNChronize!!!
 %fprintf(fgen, 'SYST:ERR?');
fprintf(fgen, 'FORM:BORD SWAP');  %configure the box to correctly accept the binary arb points
fprintf(fgen,'DATA:ARBitrary:FORMAT ABAB');

%   fprintf(fgen, 'SOURce2:COMBine:FEED: CH2'); % UNION CH1 AND CH2

%update waitbar
waitbar(.1,h,mes);

%% send I data to 33522 and setup channel 1
iBytes=num2str(length(I) *4 ); %# of bytes
headerI= ['SOURce1:DATA:ARBitrary IDATA, #' num2str(length(iBytes)) iBytes]; %create header
%header= ['SOURce1:DATA:ARBitrary IDATA, #' num2str(length(iBytes)) iBytes]; %create header
binblockBytesI = typecast(I, 'uint8');  %convert datapoints to binary before sending
% [headerI binblockBytesI]
fwrite(fgen, [headerI binblockBytesI], 'uint8'); %combine header and datapoints then send to instrument

fprintf(fgen, '*WAI');   %Make sure no other commands are executed until arb is done downloadin
%Set desired configuration for channel 1
%update waitbar
waitbar(.4,h,mes);
fprintf(fgen,'SOURce1:FUNCtion:ARBitrary IDATA'); % set current arb waveform to defined arb testrise
fprintf(fgen,'MMEM:STOR:DATA1 "INT:\IDATA.arb"');%store arb in intermal NV memory

        ON_OFF_FILTER_CH1 = ['SOURce1:FUNCtion:ARBitrary:FILTer ', (CH1_filter)];
        fprintf(fgen, ON_OFF_FILTER_CH1); % ON OFF filter


fprintf(fgen,'SOURce1:FUNCtion ARB'); % turn on arb function

command = ['SOURCE1:VOLT ' num2str(amp)]; %create amplitude command
fprintf(fgen,command); %send amplitude command

OFF1 = ['SOURCE1:VOLT:OFFSET ' num2str(offset1)]; %create amplitude command
fprintf(fgen,OFF1); % set offset to 0 V
Z1 = ['OUTPUT1:LOAD ' num2str(load1)];
fprintf(fgen,Z1); % set output load to 50 ohms

command = ['SOURCE1:FUNCtion:ARB:SRATe ' num2str(fd)]; %create sample rate command
fprintf(fgen,command);%set sample rate
%         fprintf(fgen, 'SOURce1:FUNCtion:ARBitrary: SYNChronize'); % UNION CH1 AND CH2
%Enable Output for channel 1
fprintf(fgen,CH1_source); %send ON OFF command
fprintf('"I" signal downloaded to channel 1\n\n')
%update waitbar
waitbar(.55,h,mes);

        %% -------send Q data to 33522 and setup channel 2
        qBytes=num2str(length(Q) * 4); %# of bytes
        headerQ= ['SOURce2:DATA:ARBitrary QDATA, #' num2str(length(qBytes)) qBytes]; %create header
        binblockBytesQ = typecast(Q, 'uint8');  %convert datapoints to binary before sending
        % [headerQ binblockBytesQ]
        fwrite(fgen, [headerQ binblockBytesQ], 'uint8'); %combine header and datapoints then send to instrument

        fprintf(fgen, '*WAI');   %Make sure no other commands are exectued until arb is done downloading
        %update waitbar
        waitbar(.85,h,mes);
        %Set desired configuration for channel 2
        fprintf(fgen,'SOURce2:FUNCtion:ARBitrary QDATA'); % set current arb waveform to defined arb testrise
        fprintf(fgen,'MMEM:STOR:DATA2 "INT:\QDATA.arb"'); %store arb in non V memory
        
        ON_OFF_FILTER_CH2 = ['SOURce2:FUNCtion:ARBitrary:FILTer ', (CH2_filter)];
        fprintf(fgen, ON_OFF_FILTER_CH2); % ON OFF filter
        
        fprintf(fgen,'SOURce2:FUNCtion ARB'); % turn on arb function
        command = ['SOURCE2:VOLT ' num2str(amp)]; %create amplitude command
        fprintf(fgen,command); %send amplitude command

        OFF2 = ['SOURCE2:VOLT:OFFSET ' num2str(offset2)]; %create amplitude command
        fprintf(fgen,OFF2); % set offset to 0 V
        Z2 = ['OUTPUT2:LOAD ' num2str(load2)];
        fprintf(fgen,Z2); % set output load to 50 ohms 

        command = ['SOURCE2:FUNCtion:ARB:SRATe ' num2str(fd)]; %create sample rate command
        fprintf(fgen,command);%set sample rate

         %Enable Output for channel 2
        fprintf(fgen,CH2_source); %send ON OFF command        
        fprintf('"Q" signal downloaded to channel 2\n\n')
       

%Phase sync both channels
fprintf(fgen,'PHAS:SYNC');
%get rid of message box
waitbar(1,h,mes);
delete(h);
%Read Error
fprintf(fgen, 'SYST:ERR?');
errorstr = fscanf (fgen);
% error checking
if strncmp (errorstr, '+0,"No error"',13)
   errorcheck = 'Arbitrary waveform generated without any error\n';
   fprintf (errorcheck)
else
   errorcheck = ['Error reported: ', errorstr];
   fprintf (errorcheck)
end

fclose(fgen);


     
%     % I=real(idata);
%     % I=I.';
%     % 
%     % Q=imag(qdata);
%     % Q=Q.';
%     IQData = Frame_time;
% %     figure
% %     plot(abs(fft(IQData)));
%      
%     % Define a filename for the data in the ARB
%     ArbFileName = 'MATLAB_WFM';
%      
%      
%     % Open a VISA connection or a raw TCPIP/GPIB connection to the instrument
%     % deviceObject = visa('agilent','TCPIP0::A-N5182A-80056.dhcp.mathworks.com::inst0::INSTR');
%     % deviceObject = gpib('agilent',8,19);
%     % deviceObject = tcpip('192.168.10.5',5025);
%      
%     % Find a VISA-USB object.
%     deviceObject = instrfind('Type', 'visa-usb', 'RsrcName', VISA, 'Tag', '');
%      
%     % Create the VISA-USB object if it does not exist
%     % otherwise use the object that was found.
%     if isempty(deviceObject)
%         deviceObject = visa('KEYSIGHT', VISA);
%     else
%         fclose(deviceObject);
%         deviceObject = deviceObject(1);
%     end
%      
%     % Connect to instrument object, obj1.
%      
%     % Set up the output buffer size to hold at least the number of bytes we are
%     % transferring
%     deviceObject.OutputBufferSize = 10000000;
%     % Set output to Big Endian with TCPIP objects, because we do the interleaving 
%     % and the byte ordering in code. For VISA or GPIB objecs, use littleEndian.
%     deviceObject.ByteOrder = 'bigEndian';
%      
%     % Adjust the timeout to ensure the entire waveform is downloaded before a
%     % timeout occurs
%     deviceObject.Timeout = 10.0;
%      
%     % Open connection to the instrument
%     fopen(deviceObject);
%      
%     % Seperate out the real and imaginary data in the IQ Waveform
%     wave = [real(IQData);imag(IQData)];
%     wave = wave(:)';    % transpose the waveform
%      
%     % Scale the waveform
%     tmp = max(abs([max(wave) min(wave)]));
%     if (tmp == 0)
%         tmp = 1;
%     end
%     % ARB binary range is 2's Compliment -32768 to + 32767
%     % So scale the waveform to +/- 32767 not 32768
%     scale = 2^15-1;
%     scale = scale/tmp;
%     wave = round(wave * scale);
%     modval = 2^16;
%     % Get it from double to unsigned int and let the driver take care of Big
%     % Endian to Little Endian for you  Look at ESG in Workspace.  It is a
%     % property of the VISA driver.
%     wave = uint16(mod(modval + wave, modval));
%      
%     % Some more commands to make sure we don't damage the instrument
%     fprintf(deviceObject,':OUTPut:STATe OFF')
%     fprintf(deviceObject,':SOURce:RADio:ARB:STATe ON')
%     fprintf(deviceObject,':OUTPut:MODulation:STATe ON')
%      
%     % freq =100e6;  
%     % fprintf('%d', num);
%     % Set the instrument source freq
%     fprintf(deviceObject, 'RADio:ARB:SCLock:RATE 10e6');
%     fprintf(deviceObject, 'SOURce:FREQuency 2220e6');
%     % % Set the source power
%     fprintf(deviceObject, 'POWer -25');
%      
%     % Write the data to the instrument
%     n = size(wave);
%     sprintf('Starting Download of %d Points\n',n(2)/2)
%     binblockwrite(deviceObject,wave,'uint16',[':MEM:DATA:UNProtected "WFM1:' ArbFileName '",']);
%     % Write out the ASCII LF character
%     fprintf(deviceObject,'');
%      
%     % Wait for instrument to complete download
%     % If you see a "Warning: A timeout occurred before the Terminator was reached." 
%     % warning you will need to adjust the deviceObject.Timeout value until no
%     % warning results on execution
%     commandCompleted = query(deviceObject,'*OPC?');
%      
%     % Some more commands to start playing back the signal on the instrument
%     fprintf(deviceObject,':SOURce:RADio:ARB:STATe ON')
%     fprintf(deviceObject,':OUTPut:MODulation:STATe ON')
%     fprintf(deviceObject,':OUTPut:STATe ON')
%     fprintf(deviceObject,[':SOURce:RADio:ARB:WAV "ARBI:' ArbFileName '"']);
%      
%     % Close the connection to the instrument
%     fclose(deviceObject); delete(deviceObject); clear deviceObject

end