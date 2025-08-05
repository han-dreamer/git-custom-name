function [Score, Rank, Z_normalized] = TOPSIS(X, W, K)
% TOPSIS Method for Multiple-Criteria Decision Making
%
% This function calculates the TOPSIS score and rank for a set of alternatives
% based on multiple criteria.
%
% Syntax: [Score, Rank, Z_normalized] = TOPSIS(X, W, K)
%
% Inputs:
%   X - The decision matrix (n x m), where n is the number of alternatives
%       and m is the number of criteria.
%   W - The weight vector for criteria (1 x m). Sum of weights should be 1.
%   K - The cost/benefit indicator vector (1 x m). 
%       Use 1 for benefit criteria (the larger, the better).
%       Use -1 for cost criteria (the smaller, the better).
%
% Outputs:
%   Score - The final TOPSIS score for each alternative (n x 1).
%   Rank - The rank of each alternative based on the score (n x 1).
%   Z_normalized - The normalized decision matrix after all processing.
%
% Example:
%   X = [0.1 0.2 0.3 0.4; 0.5 0.6 0.7 0.8; 0.9 1.0 1.1 1.2]; % 3 alternatives, 4 criteria
%   W = [0.25 0.25 0.25 0.25]; % Equal weights
%   K = [1 -1 1 1]; % Criteria 2 is a cost criterion
%   [S, R] = TOPSIS(X, W, K);
%
% Author: Your Name
% Date: 2024-XX-XX

%% Step 0: Input Validation
[n, m] = size(X);
if m ~= length(W) || m ~= length(K)
    error('The number of columns in X, W, and K must be the same.');
end

%% Step 1: Forward Normalization (Cost to Benefit)
% Make all criteria benefit-type
disp('Step 1: Forward Normalization...');
X_forward = X;
for j = 1:m
    if K(j) == -1
        % For cost criteria, use reciprocal or max-value subtraction
        % Here we use max-value subtraction for simplicity
        X_forward(:, j) = max(X(:, j)) - X(:, j);
    end
end

%% Step 2: Vector Normalization
disp('Step 2: Vector Normalization...');
Z_normalized = X_forward ./ sqrt(sum(X_forward.^2));

%% Step 3: Weighted Normalization
disp('Step 3: Weighted Normalization...');
V = Z_normalized .* W;

%% Step 4: Determine Ideal and Negative-Ideal Solutions
disp('Step 4: Determine Ideal and Negative-Ideal Solutions...');
Z_plus = max(V, [], 1);  % Ideal solution (Z+)
Z_minus = min(V, [], 1); % Negative-ideal solution (Z-)

%% Step 5: Calculate Distances
disp('Step 5: Calculate Distances...');
D_plus = sqrt(sum((V - Z_plus).^2, 2));   % Distance to ideal solution
D_minus = sqrt(sum((V - Z_minus).^2, 2)); % Distance to negative-ideal solution

%% Step 6: Calculate the Final Score
disp('Step 6: Calculate the Final Score...');
Score = D_minus ./ (D_plus + D_minus);

% Handle division by zero case if D_plus + D_minus is zero
if any(isnan(Score))
    disp('Warning: Some scores resulted in NaN, likely due to D+ and D- both being zero. Setting these scores to 0.');
    Score(isnan(Score)) = 0;
end

%% Step 7: Rank the Alternatives
disp('Step 7: Ranking...');
[~, Rank] = sort(Score, 'descend');

disp('TOPSIS calculation completed.');

end