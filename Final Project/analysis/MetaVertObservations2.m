function MetaVertObservations2()
    % Load datasets
    data = readtable('midtermdata-dev.csv (LastWeek_prediction).csv');
    lastWeek = string(data{:,6}) == 'TRUE'; % Whether given row is a last week
    predLastWeek = string(data{:,7}) == 'TRUE';
    
    incorrect = data{predLastWeek & ~lastWeek, 8}; % Distribution of weeks where lastWeek was pred TRUE but actually FALSE
    correct = data{predLastWeek & lastWeek, 8}; % Distribution of weeks where lastWeek was pred correctly

    figure();
    histogram(incorrect);
    xlabel('Standardized Number of Posts');
    ylabel('Number of Instances');
    title("Distribution of NumPosts where LastWeek was predicted TRUE but actually FALSE");
    
    figure();
    histogram(correct);
    xlabel('Standardized Number of Posts');
    ylabel('Number of Instances');
    title("Distribution of NumPosts where LastWeek was predicted as TRUE and actually was TRUE");
end