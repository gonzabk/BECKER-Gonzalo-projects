function [mat_compY, mat_compU, mat_compV, frameCount ] = extractFrames( fid )
%extractFrames extracts components Y, U and V of each frame in fid into
%mat_compY, mat_compU and mat_compV and returns the frame count.
%   Detailed explanation goes here

mat_compY = [] ;
mat_compU = [] ;
mat_compV = [] ;

frameCount = 0 ;

while 1
    [compY,compU,compV]=yuv_readimage(fid);
    
    if isempty(compY)
        break
    end
    
    frameCount = frameCount +1 ;
    
    mat_compY(:,:,frameCount) = compY ;
    mat_compU(:,:,frameCount) = compU ;
    mat_compV(:,:,frameCount) = compV ;

end

