clear 
close all

% Filename
file = 'news.qcif';
% Quality coefficient
Q = 50;

% Open the file
fid = fopen(file,'r');

% Extract frames from the video
[matY, matU, matV, fCount] = extractFrames(fid);

% Preallocate cell array to store decompressed frames
decompressed_frames = cell(1, fCount);

% Compression and Decompression Loop
for frame = 1:fCount
    % Get the current frame
    Image = matY(:,:,frame);
    
    % Upsample U and V matrices to match the size of the Y matrix
    ImageU = imresize(matU(:,:,frame), size(Image));
    ImageV = imresize(matV(:,:,frame), size(Image));

    % JPEG compression
    [coded_im, rs, rs_pairs, HK, HL, Q_mat, pixelCount] = JPEG(Image, Q);

    % JPEG decompression
    decoded_im = JPEG_Decomp(coded_im, rs, rs_pairs, HK, HL, Q_mat, pixelCount);
    
    % Add U, V, then reconvert to RGB
    YUV = cat(3, decoded_im, ImageU, ImageV);
    im_decomp = ycbcr2rgb(uint8(YUV));
    
    % Store decompressed frame in cell array
    decompressed_frames{frame} = im_decomp;
end

% Close the file
fclose(fid);

% Display the decompressed video
figure;
for frame = 1:fCount
    imshow(decompressed_frames{frame});
    title(['Decompressed Frame ', num2str(frame)]);
    pause(0.05); % Adjust the pause duration as needed
end
