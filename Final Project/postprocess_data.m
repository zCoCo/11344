% Post-Process the Data by Removing Outliers (using MAD) in the Features 
% which could Contain Major Outliers all Normalizing all Features.
function postprocess_data()
    %% Settings:
    % Columns to Prune Outliers From:
    pruneColumns = ["TESSMagnitude", "StarTemp", "StarRadius", "StarLuminosity", "PPPeriod", "TransitDepth", "TransitDuration", "TransitEdgeTime"];
    % Columns to Normalize (everything but class and TIC_ID):
    normColumns = ["TESSMagnitude", "StarTemp", "StarRadius", "StarLuminosity", "NFluxRegions", "PPPeriod","TransitDepth","TransitDuration","TransitEdgeTime", "PlanetRadius", "StarGap"];
    
    %% Setup:
    % Load transit dataset
    data = readtable('./data/data_with_metaparams.csv');
    
    %% Replace Numeric Class with Nominal Class (necessary for Weka):
    toi = data.class == 1;
    data.class = string(num2str(data.class)); % unfortunately, this outer #string is necessary
    data.class(toi)= "toi";
    data.class(~toi) = "not";
    
    %% Prune Outliers
    disp("Identifying Outliers . . .");
    % Remove Outliers (instances where any column contains an extreme outlier within its respective column vector):
    outlier = zeros(height(data),1);
    for c = pruneColumns
        col = data.(c);
        test = find(col ~= 0); % Don't test cells where the value is 0 (default value / N/a. value)
        [~, columnOutlier] = rmoutliers(col(test), 'median'); % Use MAD
        
        new_outlier = false(size(outlier));
        new_outlier(test(columnOutlier),1) = true;
        outlier = outlier | new_outlier;
    end
    
    % Have Human Make Sure this Won't Prune Too Much of the Dataset OR
    % Have a Bias for Pruning out Objects of Either Class (especially TOI),
    % if it were, this would mean that it's likely treating objects with
    % distinguishing features as outliers.
    disp("Identified " + sum(outlier) + " Outliers: ");
    N_toi = sum(data{~outlier,1} == "toi");
    N_non = length(data{~outlier,1}) - N_toi;
    disp("The remaining distribution is " + N_toi + " TOIs and " + N_non + " Non-TOIs after a TOI:Non-TOI removal ratio of: " + sum(data{outlier,1} == "toi")/sum(data{outlier,1} ~= "toi") + ".");
    
    removeOutliers = input("Remove Outliers (y/n)?", 's');
    if contains(removeOutliers, 'y') || contains(removeOutliers, 'Y')
        disp("Removing Outliers . . .");
        data = data(~outlier,:);
    end
    
    %% Normalize:
    disp("Normalizing Data . . .");
    for c = normColumns
        col = data.(c);
        min_val = min(col);
        max_val = max(col);
        data.(c) = (col-min_val) ./ (max_val-min_val);
    end
    
    %% Export Data:
    writetable(data, './data/complete_data.csv');
end