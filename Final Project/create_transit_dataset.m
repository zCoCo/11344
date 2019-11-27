% Processes photometric data to add transit data to the star_data dataset.
function create_transit_dataset()
    % Load star dataset
    in.data = readtable('./data/star_data.csv');
    % Load Transit Period Table
    transitPeriods = readtable(char("./data/transit_periods.csv"));

    % Input Data Column Indices:
    in.class = 1;
    in.TICID = 2;
    in.mag = 3;
    in.Teff = 4; % [K]
    in.Rstar = 5; % [Rsol]
    in.Lstar = 6; % [Lsol]

    N_peaks = zeros(height(in.data),1); % Number of Flux Histogram Peaks for Each Light Curve
    peakSep = zeros(height(in.data),1); % Mean Peak Separatation in Flux Histogram for Each Light Curve
    PPPeriod = zeros(height(in.data),1); % Peak Power Period of Object
    
    validObjects = false(height(in.data),1); % Boolean List of All TIC Objects for which All Analysis could be Performed (all data valid).
    timer = tic;
    
    for i = 1:height(in.data)
        TIC = in.data{i,in.TICID};

        if ~mod(i,10)
            DT = toc(timer);
            disp("Reading Light Curve " + i + " for TIC " + TIC);
            if i > 1
                ETA = (height(in.data)-i) * DT / (i-1);
                disp("Time Elapsed: " + floor(DT/60) + "min, ETA: " + ceil(ETA/60) + "min");
            end
        end

        valid = true;
        %% Light Curve Analysis
        % Grab Light Curve:
        try
            curve = readtable(char("./light_curves/lc_"+TIC+".csv"));
        catch e
            valid = false;
            disp("No Valid Light Curve at i="+i+", TIC="+TIC);
            warning(e.message);
        end
        % Analyze Light Curve:
        if valid
            % Find most common collections of fluxes:
            [N,X] = hist(curve.flux);
            [~,flux_clumps] = findpeaks(N,X);
            N_peaks(i) = length(flux_clumps);
            sep = mean(diff(flux_clumps));
            if isnan(sep)
                peakSep(i) = 0;
            else
                %% TODO: Precalc F0, (and DF?) before descent.
                peakSep(i) = 1e6 * sep / max(flux_clumps);
            end
        end
        
        %% Periodogram Results (peak power period)
        % Find Transit Period (if in table and TIC hasn't been invalidated yet):
        if valid
            period = transitPeriods(transitPeriods.TICID == TIC, :).Period{:};
            if numel(period) > 0
                PPPeriod(i) = period{1};
            else
                valid = false; % Period not in table (periodogram couldn't be built)
            end
        end
        
        %% Folded Transit Curve Analysis
        if valid
            
        end
        
        % Mark Validity:
        validObjects(i) = valid;
    end

    % Trim data from invalid objects:
    N_peaks = N_peaks(validObjects);
    peakSep = peakSep(validObjects);

    %% Create New Dataset:
    % Columns of Output Table (name of fields following structs,
    % 'class' must be first and is not present in structs):
    columns = {'class', 'TICID', 'mag', 'Teff', 'Rstar', 'Lstar', 'Npeaks', 'peakSep'};
    columnDesc = {'class', 'TIC_ID', 'TESSMagnitude', 'StarTemp', 'StarRadius', 'StarLuminosity', 'NFluxRegions', 'MeanFluxSep'}; % human-readable headers for output table
    %columns = {'class', 'TICID', 'mag', 'Teff', 'Rstar', 'Lstar', 'fit', 'depth', 'period', 'number', 'duration', 'riseTime'};
    %columnDesc = {'class', 'TIC_ID', 'TESSMagnitude', 'StarTemp', 'StarRadius', 'StarLuminosity', 'FitRSquared', 'DepthPPM', 'Period', 'Number', 'Duration', 'RiseTime'}; % human-readable headers for output table

    out.data = table;
    out.data{:,:} = in.data{validObjects,:}; % copy over valid entries
    out.data{:,end+1:end+2} = [N_peaks, peakSep];

    % Export Dataset:
    out.data.Properties.VariableNames = columnDesc; % Label Table
    writetable(out.data, './data/data_with_transits.csv');
end
