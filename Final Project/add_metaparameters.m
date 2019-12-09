% Add Metaparameters to the Dataset (determined from performing error
% analysis).
function add_metaparameters()
    % Load transit dataset
    data = readtable('./data/data_with_transits.csv');

    % Add MetaParameters:
    data.StarGap = abs(data.StarTemp - 4000); % 4000K = upper bound for Red Dwarfs
    
    % Stellar Mass (in SI units):
    Lstar = data.StarLuminosity; % Shorthand
    data.StarMass = (Lstar/1.4).^(1/3.5); % Default.
    % High Mass Stars:
    region = Lstar / 32000 > 55;
    data.StarMass(region) = Lstar(region) / 32000;
    % Main Sequence Stars:
    region = (0.43 <= Lstar.^0.25) & (Lstar.^0.25 < 2);
    data.StarMass(region) = Lstar(region).^0.25;
    % Red Dwarfs:
    region = (Lstar/0.23).^(1/2.3) < 0.43;
    data.StarMass(region) = (Lstar(region)/0.23).^(1/2.3);
    
    % Export Data:
    writetable(data, './data/data_with_metaparams.csv');
end