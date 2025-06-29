function [rs_vec, c_vec] = RLE(data_vec)
%Performs RLE [(r, s), c] on data_vec and outputs binary vector output_vec
rs_vec=[];
c_vec = [];

count_zeros = 0 ;
for i = 1:length(data_vec)
    if data_vec(i) == 0 
        count_zeros = count_zeros + 1 ;
    else
        % Encode sign bit separately
        sign_bit = data_vec(i) < 0;
        % Take absolute value for RLE encoding
        abs_value = abs(data_vec(i));
        % Convert to binary vector
        bin_vec = de2bi(abs_value, 'left-msb');
        % Add sign bit as the most significant bit
        bin_vec = [sign_bit bin_vec];
        
        rs_vec = [rs_vec count_zeros length(bin_vec)] ;
        count_zeros = 0 ;
        c_vec = [c_vec bin_vec] ;
    end
end
rs_vec = [rs_vec 0 0] ;




