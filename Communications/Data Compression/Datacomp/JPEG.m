function [coded_im, rs, rs_pairs, HK, HL, Q_mat, pixelCount] = JPEG(Image, Q, NBCOMP)
%For the moment only works on Y component (luminance)

%% compression
% Store the pixel count
pixelCount = numel(Image);

% perform DCT
im_dct = bdct(Image, [8, 8]) ;

% quantization
Quant_mat = csvread('Quant_mat.csv') ;
% Apply the quality coefficient Q
Q_mat = Quant_mat + (1 - Quant_mat) * (Q / 100);
% 396 is the number of 8x8 blocks in the image
im_quant = round(im_dct ./ (Q_mat(:) * ones(1,396))) ; % ATTENTION

selected_comp = zeros(size(im_quant)) ;
selected_comp(1:NBCOMP,:) = im_quant(1:NBCOMP,:) ;

% ZIGZAG
    % Initialize an empty vector for the result
    [rows, cols] = size(im_quant);
    im_dct_flat = [];
    
    % Iterate over each block (here each colum in im_dct represents an 8x8
    % block
    for r = 1:cols
        %Reshape the colum into an 8x8 block
        block = reshape(selected_comp(:,r),8,8);
        im_dct_flat = [im_dct_flat zigzag_scan(block)];
    end


% method 2 (from lab teacher)
[rs, c] = RLE(im_dct_flat) ;

% Huffman coding on resulting sequence
[rs_pairs, proba] = CountProbabilities(rs) ;
HL = HuffLen(proba) ;
HLlen = HuffTabLen(HL) ;
HK = HuffCode(HL, 1);

coded_im = [] ;
c_pointer = 1 ;

for idx = 1:2:length(rs)
    r = rs(idx) ;
    s = rs(idx+1) ;
    is_codeword = false ;
    j = 1 ;
    while ~is_codeword
        if isequal([r,s], rs_pairs(j,:))
            is_codeword = true ;
            coded_im = [coded_im HK(j, 1:HL(j))] ;
            coded_im = [coded_im c(c_pointer:c_pointer+s-1)] ;
            c_pointer = c_pointer + s;
        else
            j=j+1 ;
        end
    end
end

    


