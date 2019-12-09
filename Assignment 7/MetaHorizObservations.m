function MetaHorizObservations()
    % Configurable Parameters:
    evaluateColumn = 'Lstar'; % ID (see below) of column of feature to be observed
    columnLongName = 'Star Luminosity'; % Cosmetic long-form name of column to be observed (used in plots)
    fileAddr = 'complete_data-test.csv (class_prediction).csv';
    
    classCol = 13;
    predCol = 14;
    
    % Load datasets
    data = readtable(fileAddr);
    toi = string(data{:,classCol}) == 'toi'; % Whether given row is actually a TOI
    predTOI = string(data{:,predCol}) == 'toi'; % Whether given row was predicted to be a TOI
    
    % Numeric Column IDs:
    cols.Lstar = 5;
    
    % Find Index of Column being Evaluated:
    colIdx = cols.(char(evaluateColumn));
    
    incorrect = data{~predTOI & toi, colIdx}; % Distribution of feature where class was pred not-TOI but actually TOI
    correct = data{predTOI & toi | ~predTOI & ~toi, colIdx}; % Distribution of feature where class was pred correctly

    figure();
    histogram(incorrect);
    xlabel(char("Feature: " + columnLongName + " (" + evaluateColumn + ")"));
    ylabel('Number of Instances');
    title("Distribution of " + evaluateColumn + " where class was predicted as Non-TOI but actually TOI");
    
    figure();
    histogram(correct);
    xlabel(char("Feature: " + columnLongName + " (" + evaluateColumn + ")"));
    ylabel('Number of Instances');
    title("Distribution of " + evaluateColumn + " where class was predicted correctly");
end