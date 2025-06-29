clear 
close all

NBCOMP = 64 ;
Q = 100; % en %
FILE = 'news.qcif';

%% Load and display image

fprintf('\n----------------------------\n');
fprintf('File: %s\n', FILE);
fprintf('Nb of components kept: %d\n', NBCOMP);
fprintf('Q = %d %% \n', Q);
fprintf('----------------------------\n');

% Open the file
fid = fopen(FILE,'r');

[matY, matU, matV, fCount] = extractFrames(fid);

Image = matY(:,:,1);


imagesc(Image);
colormap(gray)

%% Compression

% perform DCT
im_dct = bdct(Image, [8, 8]) ;
figure 
%imagesc(im_dct) ;


% quantization
Q_mat = csvread('Quant_mat.csv') ;
Q_mat = Q_mat + (ones(size(Q_mat))-Q_mat)/100*Q;
im_quant = round(im_dct ./ (Q_mat(:) * ones(1,396))) ; % ATTENTION: not portable
                                                       % bcs hard coded img
                                                       % size
% select limited amount of compnents
selected_comp = zeros(size(im_quant)) ;
selected_comp(1:NBCOMP,:) = im_quant(1:NBCOMP,:) ;



% perform zigzag
im_dct_flat = zeros(1, numel(selected_comp));
% Zigzag of the first block
prev_zz = zigzag(reshape(selected_comp(:, 1), [8, 8]));
% Store the first block as it is
im_dct_flat(1:64) = prev_zz; 
for block = 2:size(im_dct, 2)
    curr_zz = zigzag(reshape(selected_comp(:, block), [8, 8])); 
    diff = curr_zz(1,1) - prev_zz(1,1);
    prev_zz = curr_zz;
    curr_zz(1,1)=diff;
    im_dct_flat((block - 1) * 64 + 1 : block * 64) = curr_zz; 
end

% perform RLE on im_dct_flat
[rs, c] = RLE(im_dct_flat) ;

% Huffman coding on resulting sequence
[rs_pairs, proba] = CountProbabilities(rs) ;
HL = HuffLen(proba) ;
HLlen = HuffTabLen(HL) ;
HK = HuffCode(HL);

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
im_quant_est = zeros(size(im_quant));
prev_zz = izigsc(im_dct_flat_decod(1:64), 8);
im_quant_est(:,1) = prev_zz(:);
for block = 2:size(im_dct, 2)
    curr_zz = izigsc(im_dct_flat_decod((block - 1) * 64 + 1 :64*block), 8);
    curr_zz(1,1) = curr_zz(1,1) + prev_zz(1,1);
    im_quant_est(:,block) = curr_zz(:);
    prev_zz = curr_zz;
end

% inverse DCT

im_idct = ibdct(im_quant_est,[8,8], size(Image)) ;
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
E = sum(epsilon.^2, 'all')/(M*N);
PSNR = 10*log10(N_g^2/E);
fprintf('PSNR = %f\n', PSNR);

fclose(fid);