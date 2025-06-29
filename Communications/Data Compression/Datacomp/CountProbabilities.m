function [rs_pairs, probabilities] = CountProbabilities(rs_vec)
    % Initialize containers
    rs_pairs = [];
    probabilities = [];

    % Iterate over rs_vec to count probabilities
    for i = 1:2:length(rs_vec)-1
        r = rs_vec(i);
        s = rs_vec(i+1);
        
        % Check if the pair already exists in rs_pairs
        pair_exists = false;
        for j = 1:size(rs_pairs, 1)
            if all(rs_pairs(j, :) == [r, s])
                % Increment the probability if the pair exists
                probabilities(j) = probabilities(j) + 1;
                pair_exists = true;
                break;
            end
        end
        
        if ~pair_exists
            % Add the new pair to rs_pairs and set its probability to 1
            rs_pairs = [rs_pairs; [r, s]];
            probabilities = [probabilities; 1];
        end
    end

    % Convert counts to probabilities
    total_elements = length(rs_vec) / 2; % Since rs_vec contains pairs
    probabilities = probabilities / total_elements;
end
