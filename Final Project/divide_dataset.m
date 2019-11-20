% Splits Complete Dataset into training (development), testing, and final
% testing datasets.
function divide_dataset()
    % Dataset Partitioning:
    p_train = 1/2; % fraction of complete dataset reserved for training
    p_test = 2/6; % fraction of complete dataset reserved for testing
    p_fin = 1-p_train-p_test; % fraction of complete dataset reserved for final testing

    % Load dataset:
    in = readtable('./data/complete_data.csv');
    
    % Partition Data by Index (randomly):
    h = height(in);
    idxPerm = randperm(h); % Permutation of all indices
    idx_train = idxPerm(1:floor(p_train*h)); % Training indices
    idx_test = idxPerm(floor(p_train*h)+1:floor(p_train*h)+1+floor(p_test*h)); % Testing indices
    idx_fin = idxPerm(floor(p_train*h)+floor(p_test*h)+2:end); % Final testing indices
    
    % Split up Datasets:
    train = in(idx_train,:);
    test = in(idx_test,:);
    fin = in(idx_fin,:);
    
    % Save Partitioned Tests:
    writetable(train, './data/complete_data-train.csv');
    writetable(test, './data/complete_data-test.csv');
    writetable(fin, './data/complete_data-fin.csv');
end