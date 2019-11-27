function MetaHorizObservations()
    % Configurable Parameters:
    evaluateColumn = 'mag'; % ID (see below) of column of feature to be observed
    columnLongName = 'Magnitude'; % Cosmetic long-form name of column to be observed (used in plots)

    % Load datasets
    data = readtable('midtermdata-dev.csv (LastWeek_prediction).csv');
    toi = string(data{:,6}) == 1; % Whether given row is actually a TOI
    predTOI = string(data{:,7}) == 1; % Whether given row was predicted to be a TOI
    
    % Numeric Column IDs:
    cols.TICID = 2;
    cols.mag = 3;
    cols.Teff = 4;
    cols.Rstar = 5;
    cols.Lstar = 6;
    cols.Npeaks = 7;
    cols.peakSep = 8;
    
    % Find Index of Column being Evaluated:
    colIdx = cols.(char(evaluateColumn));
    
    incorrect = data{predTOI & ~toi, colIdx}; % Distribution of feature where class was pred TOI but actually NonTOI
    correct = data{predTOI & toi | ~predTOI & ~toi, colIdx}; % Distribution of feature where class was pred correctly

    figure();
    histogram(incorrect);
    xlabel(char("Feature: " + columnLongName + " (" + evaluateColumn + ")"));
    ylabel('Number of Instances');
    title("Distribution of " + evaluateColumn + " where class was predicted as TOI but actually Non-TOI");
    
    figure();
    histogram(correct);
    xlabel(char("Feature: " + columnLongName + " (" + evaluateColumn + ")"));
    ylabel('Number of Instances');
    title("Distribution of " + evaluateColumn + " where class was predicted correctly");
end