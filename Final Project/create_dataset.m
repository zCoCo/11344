% Used to create the full dataset of TOIs and non-TOIs by combining the
% list of ~1k TOIs with twice as many entries from the CTL with the 
% highest priority that are not in the list of TOIs.
function create_dataset()
    % Load datasets
    toi.data = readtable('./raw_data/TOI.csv');
    ctl.data = readtable('./raw_data/ObservedCTL.csv');

    % Columns of Output Table (name of fields following structs, 
    % 'class' must be first and is not present in structs):
    columns = {'class', 'ticID', 'mag', 'Teff', 'Rstar', 'Lstar'};
    columnDesc = {'class', 'TIC_ID', 'TESSMagnitude', 'StarTemp', 'StarRadius', 'StarLuminosity'}; % human-readable headers for output table
    
    % TOI Column Indices:
    toi.ticID = 1;
    toi.toiID = 2;
    toi.mag = 14;
    toi.RA = 19;
    toi.dec = 20;
    toi.Teff = 42; % [K]
    toi.logg = 44; % [cm/s^2]
    toi.Rstar = 46; % [Rsol]
    toi.sector = 48;
    % CTL Column Indices:
    ctl.RA = 1;
    ctl.dec = 2;
    ctl.mag = 3;
    ctl.Teff = 4; % [K]
    ctl.Rstar = 6; % [Rsol]
    ctl.Mstar = 7; % [Msol]
    ctl.sector = 10;
    ctl.ticID = 27;
    ctl.Lstar = 31; % [Lsol]
    
    % Fill in missing but solveable toi features:
    toi.Lstar = width(toi.data)+1;
    Tsol = 5578; % Surface temperature of sun [K]
    toi.data{:,toi.Lstar} = toi.data{:,toi.Rstar}.^2 .* (toi.data{:,toi.Teff}./Tsol).^4; % Solar Stefan-Boltzman Law
    
    % Build indices of columns to be selected from CTL table in order of
    % columns vector:
    toi.columns = zeros(numel(columns)-1,1);
    for j = 2:numel(columns)
        c = columns{j};
        toi.columns(j-1) = toi.(c);
    end
    
    % Build indices of columns to be selected from CTL table in order of
    % columns vector:
    ctl.columns = zeros(numel(columns)-1,1);
    for j = 2:numel(columns)
        c = columns{j};
        ctl.columns(j-1) = ctl.(c);
    end
    
    % Tells whether the data in struct s at row r and the column given by 
    % field c in struct s is valid (this centralizes the definition of 
    % valid).
    function v = isValid(s, r,c)
        v = ~isnan(s.data{r,s.(c)}) & s.data{r,s.(c)} > 0;
    end
    
    % Collect Indices of Valid TOI entries:
    validTOI = []; % indices of valid TOI entries
    i = 1;
    while i <= height(toi.data)
        % Ensure all required columns are not NaN and >0
        valid = true;
        for j = 2:numel(columns)
            c = columns{j};
            valid = valid & isValid(toi, i,c);
        end
        if valid
            validTOI(end+1) = i;
        end
        i = i+1;
    end
    
    % Collect Indices of Valid non-TOI CTL entries:
    validCTL = []; % indices of valid CTL entries
    i = 1;
    while numel(validCTL) < 2 * height(toi.data) && i <= height(ctl.data)
        % Ensure all required columns are not NaN and >0
        valid = true;
        for j = 2:numel(columns)
            c = columns{j};
            valid = valid & isValid(ctl, i,c);
        end
        if valid
            validCTL(end+1) = i;
        end
        i = i+1;
    end
    
    % Start table with all valid TOI entries:
    n = numel(validTOI);
    out = table();
    out{1:n,:} = [ones(n,1), toi.data{validTOI, toi.columns}];
    
    % Append Valid CTL Entries to Table:
    n = numel(validCTL);
    out{end+1:end+n,:} = [zeros(n,1), ctl.data{validCTL, ctl.columns}];
    
    % Shuffle Dataset:
    out{:,:} = out{randperm(height(out)), :};
    
    % Export Dataset:
    out.Properties.VariableNames = columnDesc; % Label Table
    writetable(out, './data/star_data.csv');
end