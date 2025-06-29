% This function implements an inverse zigzag scan of a 1D array as an input.
%
% Syntax:  [A] = izigsc(B,dim)
%
% Inputs:
%    B - Input 1D data array. The length of B should be equal to dim*dim.
%    dim - Dimension of the square zigzagged matrix.
%
% Output:
%    A - Outcome matrix after inverse zigzag scanning.
% 
% Example: 
%    B = [1:16];
%    dim = 4;
%    A = [1,2,6,7;3,5,8,13;4,9,12,14;10,11,15,16];
%
%
% Author: Maryam Sadeghi
% May 2020.
%-------------------- BEGIN OF CODE ---------------------
function [A] = izigsc(B,dim)
v = ones(1,dim); k = 1;
A = zeros(dim,dim);
for i = 1:2*dim-1
    C1 = diag(v,dim-i);
    C2 = flip(C1(1:dim,1:dim),2);
    C3 = B(k:k+sum(C2(:))-1);
    k = k + sum(C2(:));
    if mod(i,2) == 0
       C3 = flip(C3);
    end
    C4 = zeros(1,dim-size(C3,2));
    if i >= dim
       C5 = cat(2,C4, C3); 
    else       
        C5 = cat(2,C3,C4);
    end
    C6 = C2*diag(C5);
    A = C6 + A;
end