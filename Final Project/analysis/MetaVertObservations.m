function MetaVertObservations()
    % Load datasets
    data = readtable('midtermdata-dev.csv (LastWeek_prediction).csv');
    lastWeek = string(data{:,6}) == 'TRUE'; % Whether given row is a last week
    predLastWeek = string(data{:,7}) == 'TRUE';
    
    incorrect = data{predLastWeek & ~lastWeek, 5}; % Distribution of weeks where lastWeek was pred TRUE but actually FALSE
    correct = data{predLastWeek & lastWeek, 5}; % Distribution of weeks where lastWeek was pred correctly

    figure();
    histogram(incorrect);
    xlabel('Standardized In-Degree');
    ylabel('Number of Instances');
    title("Distribution of in-degree where LastWeek was predicted TRUE but actually FALSE");
    
    figure();
    histogram(correct);
    xlabel('Standardized In-Degree');
    ylabel('Number of Instances');
    title("Distribution of in-degree where LastWeek was predicted as TRUE and actually was TRUE");
end