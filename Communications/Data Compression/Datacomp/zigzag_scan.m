function zigzag_vector = zigzag_scan(matrix)
    % Get the size of the input matrix
    [rows, cols] = size(matrix);
    
    % Initialize the zigzag vector
    zigzag_vector = zeros(1, rows*cols);
    
    % Initialize indices for traversing the zigzag pattern
    idx = 1;
    row = 1;
    col = 1;
    
    % Flag to indicate the direction of traversal
    going_up = true;
    
    % Zigzag scan the matrix
    while idx <= numel(matrix)
        % Store the current element in the zigzag vector
        zigzag_vector(idx) = matrix(row, col);
        
        % Move to the next element based on the direction of traversal
        if going_up
            if col == cols
                % If at the last column, move down
                row = row + 1;
                going_up = false;
            elseif row == 1
                % If at the first row, move right
                col = col + 1;
                going_up = false;
            else
                % Otherwise, move diagonally up-right
                row = row - 1;
                col = col + 1;
            end
        else
            if row == rows
                % If at the last row, move right
                col = col + 1;
                going_up = true;
            elseif col == 1
                % If at the first column, move down
                row = row + 1;
                going_up = true;
            else
                % Otherwise, move diagonally down-left
                row = row + 1;
                col = col - 1;
            end
        end
        
        % Increment the index
        idx = idx + 1;
    end
end
