% Add Metaparameters to the Dataset (determined from performing error
% analysis).
function add_metaparameters()
    % Load transit dataset
    data = readtable('./data/data_with_transits.csv');

    % Add MetaParameters:
    data.StarGap = abs(data.StarTemp - 4000); % 4000K = upper bound for Red Dwarfs
    
    % Export Data:
    writetable(data, './data/data_with_metaparams.csv');
end