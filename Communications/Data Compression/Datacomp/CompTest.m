clear 
close all

% Filename
file = 'news.qcif';
%Quality coefficient
Q = 50;

% Open the file
fid = fopen(file,'r');

[matY, matU, matV, fCount] = extractFrames(fid);

Image = matY(:,:,1);
% Upsample U and V matrices to match the size of the Y matrix
ImageU = imresize(matU(:,:,1), size(Image));
ImageV = imresize(matV(:,:,1), size(Image));

% Stack Y, U, and V channels to form the YUV image
YUV_image = cat(3, Image, ImageU, ImageV);

% Convert YUV image to RGB color space
RGB_image = ycbcr2rgb(uint8(YUV_image));  % Convert to uint8 for correct conversion

% Display the original RGB image
figure
imagesc(RGB_image);
title('Original Image');

%JPEG compression
[coded_im, rs, rs_pairs, HK, HL, Q_mat, pixelCount] = JPEG(Image, Q);

%JPEG decompression
decoded_im = JPEG_Decomp(coded_im, rs, rs_pairs, HK, HL, Q_mat, pixelCount);

% Add U, V then reconvert to RGB
YUV = cat(3, decoded_im, ImageU, ImageV);
RGB_Decomp = ycbcr2rgb(uint8(YUV));

%Display
figure
imagesc(RGB_Decomp);
title('Compressed Image');