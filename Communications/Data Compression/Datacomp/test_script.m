clear 
close all

NBCOMP = 32 ;
Q = 90;

% Filename
file = 'news.qcif';

% Open the file
fid = fopen(file,'r');

[matY, matU, matV, fCount] = extractFrames(fid);

%% playing with colors

Image = matY(:,:,1);
ImageU = matU(:,:,1);
ImageV = matV(:,:,1);

figure
imagesc(Image);
colormap(gray)
% figure
% imagesc(ImageU.');
% figure
% imagesc(ImageV.');
% 
% imUx2 = kron(ImageU,ones(2));
% imVx2 = kron(ImageV,ones(2));
% figure
% imagesc(Image - imUx2.' - imVx2.');


%% compression

% perform DCT
im_dct = bdct(Image, [8, 8]) ;
figure 
%imagesc(im_dct) ;


% quantization
Q_mat = csvread('Quant_mat.csv') ;
Q_mat = Q_mat + (ones(size(Q_mat))-Q_mat)/100*Q;

im_quant = round(im_dct ./ (Q_mat(:) * ones(1,396))) ; % ATTENTION
%im_quant = round(im_dct ./ (reshape(Q_mat.',[],1) * ones(1,396))) ;
%im_quant = im_dct ;

selected_comp = zeros(size(im_quant)) ;
selected_comp(1:NBCOMP,:) = im_quant(1:NBCOMP,:) ;
%selected_comp(:,1:NBCOMP) = im_quant(:,1:NBCOMP) ;



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

% perform RLE on im_dct_flat (2 methods)

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

%% decompression

% Huffman decoding 
   

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
im_dct_flat_decod = zeros(size(im_dct_flat));
idx=1;
for i = 1:length(c_decoded)
    idx = idx + rs_pairs_decoded(i,1);
    im_dct_flat_decod(idx) = c_decoded(i);
    idx = idx + 1;
end

% reverse zigzag
% initialize the un-zigzagged matrix
zigzagnt = zeros(size(im_dct));
    for r = 1:396
        block = izigsc(im_dct_flat_decod((r-1)*64 + 1 : r*64), 8);
        zigzagnt(:,r)= block(:);
    end

% reverse quantization
quantiznt = zigzagnt .* (Q_mat(:) * ones(1,396));

% IDCT
im_idct = ibdct(quantiznt,[8,8], size(Image)) ; 

imagesc(im_idct)
colormap(gray)

%% metrics

fprintf('METRICS:\n');

% compression rate (volume origibnal / compressed)
[M, N] = size(Image) ;
nb_bits_ppix_comp = numel(coded_im)/(M*N);
tx_comp = 8/nb_bits_ppix_comp;
fprintf('Tx compression = %f\n', tx_comp);

% bir rate
fprintf('Debit = %f\n', nb_bits_ppix_comp);

% PSNR
N_g = max(Image(:));
epsilon = Image - im_idct;
E = sum(epsilon(:).^2)/(M*N);
PSNR = 10*log10(N_g^2/E);
fprintf('PSNR = %f\n', PSNR);

fclose(fid);