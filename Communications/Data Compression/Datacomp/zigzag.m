function C = zigzag(A)
% example : zigzag(reshape((1:100), [10,10]).')
r = size(A,1); % number of rows of the input matrix
c = size(A,2); % number of columns of the input matrix
kk = 2;
C = [];
while kk <= r+c % the lowermost diagonal has only one element
                % which has the index r by c
                % scan the matrix as long as the last
                % diagonal is retrieved
    B = [];
    % iterate through every element of the input matrix
for ii = 1:r
    for jj = 1:c
        if ii + jj == kk % sum of the indices in a particular diagonal
                         % are the same
            B = [B,A(ii,jj)];
        end
    end
end
if mod(kk,2) == 0
    C = [C,flip(B)]; % reverse the order of the diagonal entries
                     % evenly
else
    C = [C,B];
end
kk = kk+1;
end
        
end