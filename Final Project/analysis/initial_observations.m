function initial_observations()
    % Load datasets
    raw_data = readtable('../data/complete_data-train.csv');
    toi = string(raw_data{:,1}) == "toi"; % Whether given instance is a toi
    N_toi = sum(toi) % Number of TOI Objects
    N_nontoi = height(raw_data) - N_toi % Number of Non-TOI Objects

    % Numeric Column IDs:
    cols.TICID = 2;
    cols.mag = 3;
    cols.Teff = 4;
    cols.Rstar = 5;
    cols.Lstar = 6;
    cols.Npeaks = 7;
    cols.peakSep = 8;
    
    %% Compare All Numerical Columns:
    columns = string(fields(cols));
    for j = 1:numel(columns)
        column = columns(j); % in every (numeric) column
        
        N_bins = ceil(1 + log2(min([N_toi, N_nontoi]))); % Sturges rule
        
        figure();
        [nontoi_counts, nontoi_values] = hist(raw_data{~toi, cols.(column)}, N_bins);
        [toi_counts, toi_values] = hist(raw_data{toi, cols.(column)}, N_bins);
        bar(nontoi_values, 100*nontoi_counts/N_nontoi, 'b');
        hold on
        bar(toi_values, 100*toi_counts/N_toi, 'g');
        xlabel("Feature Value: " + column, 'Interpreter','latex');
        legend({'Non-TOI Objects', 'TOI Objects'}, 'Interpreter','latex');
        ylabel('Percentage of Objects in Class', 'Interpreter','latex');
        title({"Distributions of " + column + " for TOI and Non-TOI Objects"}, 'Interpreter','latex');
        
        % Display Aggregate Results for Column:
        nontoiMean = mean(raw_data{~toi, cols.(column)});
        toiMean = mean(raw_data{toi, cols.(column)});
        fprintf("The avg. %s for Non-TOI objects was %0.3f but %0.3f for TOI objects (percent difference: %0.1f%%)\n", column, nontoiMean,toiMean, percentDiff(nontoiMean,toiMean));
    end
    
    %% Other Interesting Analysis:
    
    % Display Handful of Entries from Each:
    disp("Sample Non-TOI Objects:");
    nontoiData = raw_data(~toi, :);
    nontoiData(:,:) = nontoiData(randperm(height(nontoiData)),:); % shuffle rows
    head(nontoiData)
    disp("Sample TOI Objects:");
    toiData = raw_data(toi, :);
    toiData(:,:) = toiData(randperm(height(toiData)),:); % shuffle rows
    head(toiData)
end

function pd = percentDiff(a,b)
    pd = 100*abs(a-b) / mean([a,b]);
end