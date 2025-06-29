function im_idct = JPEG_Decomp(coded_im, rs, rs_pairs, HK, HL, Q_mat, pixelCount)



%% decompression

% Huffman Decoding (Reverse Step)
rs_pairs_decoded = zeros(length(rs)/2, 2);
c_decoded = zeros(1,length(rs)/2);
idx = 1;
i = 0;
finished=false;
while finished == false
    i = i + 1;
    for j = 1:length(HK)
        if isequal(coded_im(idx:idx+HL(j)-1), HK(j,1:HL(j)))
            rs_pairs_decoded(i,:) = rs_pairs(j,:);
            if isequal(rs_pairs_decoded(i,:), [0,0])
                finished = true;
                break;
            end
            comp_bin = coded_im(idx + HL(j) : idx + HL(j) + rs_pairs_decoded(i,2) - 1);
            c_decoded(i) = -1*(2*comp_bin(1)-1)*bi2de(comp_bin(2:end),'left-msb');
            idx = idx + HL(j) + rs_pairs_decoded(i,2);
            break;
        end
    end
end


% inverse RLE
% im_dct_flat has the same number of elements as the original image
im_dct_flat_decod = zeros(pixelCount);
idx=1;
for i = 1:length(c_decoded)
    idx = idx + rs_pairs_decoded(i,1);
    im_dct_flat_decod(idx) = c_decoded(i);
    idx = idx + 1;
end


% reverse zigzag
% initialize the un-zigzagged matrix
zigzagnt = zeros(64,396); % size(im_dct)
    for r = 1:396
        block = izigsc(im_dct_flat_decod((r-1)*64 + 1 : r*64), 8);
        zigzagnt(:,r)= block(:);
    end

% reverse quantization
quantiznt = zigzagnt .* (Q_mat(:) * ones(1,396));

% IDCT
im_idct = ibdct(quantiznt,[8,8], [144,176]) ; %size(Image)


