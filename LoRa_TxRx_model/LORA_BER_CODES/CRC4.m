function [CRC_Bits] = CRC4(input);


%% Тип CRC
% if(crc_type==4)
    crc_poly = [1, 0, 0, 1, 1];    
    crc_len  = 4;
% elseif(crc_type==3)
%     crc_poly = [1, 0, 1, 1];    
%     crc_len  = 3;
% end
crc_rem   = zeros(1,crc_len+1);
tmp_array = vertcat(input, zeros(crc_len,1));

%% Вычисление CRC
for n=1:length(input)+crc_len
    for m=1:crc_len 
        crc_rem(m) = crc_rem(m+1); 
    end;
    crc_rem(crc_len+1) = tmp_array(n);

    if(crc_rem(1) ~= 0)
        for(m=1:crc_len+1) 
            crc_rem(m) = mod(crc_rem(m)+crc_poly(m), 2); 
        end; 
    end; 
end;

    for n=1:crc_len 
        CRC_Bits(n,1) = crc_rem(n+1); 
    end;

end