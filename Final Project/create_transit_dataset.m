% Processes photometric data to add transit data to the star_data dataset.
function create_transit_dataset()
    dbstop if error; % Create debugger breakpoint so data isn't lost if there's an error at a bad time.
    % Load star dataset
    in.data = readtable('./data/star_data.csv');
    % Load Transit Period Table
    transitPeriods = readtable(char("./data/transit_periods.csv"));

    % Constants:
    Rsol = 695510e3; % [m] Radius of the Sun
    REarth = 6371e3; % [m] Radius of the Earth
    
    % Input Data Column Indices:
    in.class = 1;
    in.TICID = 2;
    in.mag = 3;
    in.Teff = 4; % [K]
    in.Rstar = 5; % [Rsol]
    in.Lstar = 6; % [Lsol]

    % Light Curve Features:
    N_peaks = zeros(height(in.data),1); % Number of Flux Histogram Peaks for Each Light Curve
    peakSep = zeros(height(in.data),1); % Mean Peak Separatation in Flux Histogram for Each Light Curve
    
    % Transit Features:
    PPPeriod = zeros(height(in.data),1); % Peak Power Period of Object
    transitDepth = zeros(height(in.data),1); % [PPM] Max Depth of Transit Event
    transitDuration = zeros(height(in.data),1); % [days] Full Duration of Transit Event
    transitEdgeTime = zeros(height(in.data),1); % [days] Time planet spends crossing the edge of the star. Average of Rise Time and Fall Time of Transit Event.
    planetRadius = zeros(height(in.data),1); % [REarth] Predicted Planet Radius
    impactParam = zeros(height(in.data),1); % Describes inclination of planet orbit w.r.t. the view from Earth.
    linearLimbDarkening = zeros(height(in.data),1); % Linear Limb-Darkening Coefficient
    quadLimbDarkening = zeros(height(in.data),1); % Quadratic Limb-Darkening Coefficient
    transitCurveFit = zeros(height(in.data),1); % Adjusted Coefficient of Determination for Fitted Transit Curve

    validObjects = false(height(in.data),1); % Boolean List of All TIC Objects for which All Analysis could be Performed (all data valid).
    timer = tic;
    
    exit_fig = figure();
    text(mean(xlim),mean(ylim), "Click this window and hold 'r' to terminate early");
    for i = 1:height(in.data)
        kkey = get(exit_fig,'CurrentCharacter');
        if strcmp(kkey, 'r')
            break;
        end
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
            disp("No Valid Light Curve at i="+i+" for TIC="+TIC);
            warning(e.message);
        end
        % Analyze Light Curve:
        if valid
            % Find most common collections of fluxes:
            [N,flux_bins] = hist(curve.flux);
            [~,flux_clumps] = findpeaks(N,flux_bins);
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
            period = transitPeriods(transitPeriods.TICID == TIC, :).Period;
            if numel(period) > 0
                period = string(period{1});
                PPPeriod(i) = sscanf(period, '%f');
            else
                valid = false; % Period not in table (periodogram couldn't be built)
            end
        end
        
        %% Folded Transit Curve Analysis
        % Grab Folded Transit Curve:
        if valid
            try
                transit = readtable(char("./transits/tr_"+TIC+".csv"));
            catch e
                valid = false;
                disp("No Valid Folded Transit Curve at i="+i+" for TIC="+TIC);
                warning(e.message);
            end
        end
        
        if valid
            % Compute Depth:
            % Bin Fluxes (TODO: Remove this once binsize in photometry collector is fixed).
            binned_flux = smooth(transit.flux, floor(numel(transit.flux)/25), 'moving');
            binned_time = transit.time; % Eventually, there could be a dimensionality reduction associated with smoothing which will make wrapping this necessary.

            % Find most common collections of fluxes in folded transit curve:
            [N,flux_bins] = hist(binned_flux);
            [~,flux_clumps] = findpeaks(N,flux_bins);
            % Include edge "peaks" if applicable (particularly the max, since this is likely the
            % baseline):
            if N(end) > N(end-1)
                flux_clumps(end+1) = flux_bins(end);
            end
            if N(1) > N(2)
                flux_clumps(end+1) = flux_bins(1);
            end

            % Nominal (baseline) Flux:
            tr_base = max(flux_clumps);
            % Maximum Transit Depth:
            tr_min = min(flux_clumps);

            % Transit Depth:
             depthRatio = (tr_base - tr_min) / tr_base;

            if depthRatio ~= 0 % Significant transit observed
             % Compute Timing:
                % Find time of first and last point below 50% depth:
                idx_50 = find(binned_flux < (tr_base + 0.5*(tr_min-tr_base))); % Index of all points below 50% transit depth
                t_50drop = transit.time(idx_50(1)); % Time when first transit drops below 50%
                t_50rise = transit.time(idx_50(end)); % Time when transit rises above 50% for last time
                % Find time of first point beyond 100% depth:
                idx_100 = find(binned_flux < (tr_base + 1.0*(tr_min-tr_base))); % Index of all points below 100% transit depth
                t_100drop = transit.time(idx_100(1)); % Time when first transit drops below 100%
                t_100rise = transit.time(idx_100(end)); % Time when transit rises above 100% for last time
                % Compute 50%->100% "Edge" Time:
                t_fall = t_100drop - t_50drop; 
                t_rise = t_50rise - t_100rise;
                transitEdgeTime(i) = 2*mean([t_fall, t_rise]) * PPPeriod(i);
                transitDuration(i) = ( (t_50rise - t_50drop) + transitEdgeTime(i) ) * PPPeriod(i);
                [~, min_idx] = min(binned_flux);
                t0 = binned_time(min_idx);
            end
        end
        
        
        % Add Features from Transit Curve-Fitting
        if valid && depthRatio ~= 0 % Significant transit observed
            % Get initial transit curve parameters:
            period = PPPeriod(i); % shorthand for later
            impactParam(i) = 0.01; % initial default for search
            transitTime = transitDuration(i);
            t0 = t0*period;
            Q0 = [depthRatio, impactParam(i), transitTime, t0, 0.5,0.25]';
            
            %disp("Fitting Transit Curve . . . ");
            % Get baseline goodness-of-fit:
            goodness_of_fit = AR2(Q0, binned_time,binned_flux, period);
            %disp(" > Pre-fit Adjusted R2: ");
            %disp(goodness_of_fit);

            if(goodness_of_fit < -10)
                % Don't even try to optimize if initial fit is absolutely horrible (to save time).
                %warning(" > Not attempting curve-fit optimization. AR2 initial is " + goodness_of_fit);
                Q = Q0;
            else
                try
                    %Q0_reduced = Q0([1 2 3 5 6]);
                    options = optimoptions('fmincon', 'Display','iter', 'Algorithm','active-set', 'MaxIterations', 30, 'TolFun', 5e-6);
                    Q = fmincon(...
                            @(qq) fit_error(qq, binned_time,binned_flux, period), ... % minimize least-squared error
                            Q0, ... % initial condition (from analyzing folded curve)
                            [0,0,0,0,1,1; ... % gam1 + gam2 < 1 -> Q5 + Q6 < 1
                             ...%0,0,1/period,2/period,0,0; ... % (t0-ttrans/2)/period >= -0.5 (not in reduced problem)
                            ],...%0,0,1/period,-2/period,0,0], ... % (t0+ttrans/2)/period <= 0.5 (not in reduced problem)
                            [0.999],...%; 1; 1], ...
                            [],[], ... % <- ignore lin. eq. constraints,
                            [0.75*depthRatio,  0, 0.75*transitTime, max(-0.5*period+transitTime/2, t0-0.25*period), 0,0], ... % lower bound (treat as tuning, keep somewhat near recovered params)
                            [1.25*depthRatio, Inf, 1.25*transitTime, min(0.5*period-transitTime/2, t0+0.25*period), 1,1], ... % upper bound (treat as tuning, keep somewhat near recovered params).
                            @Q_con, ...                                        % upper bound on impactParam is handled here
                            options ...
                        );
                    % Recreate full-Q: 
                    %Q = [Q(1); Q(2); Q(3); t0; Q(4); Q(5)];
                    goodness_of_fit = AR2(Q, binned_time,binned_flux, period);
                    disp(" > Post-fit Adjusted R2: ");
                    disp(goodness_of_fit);
                    % Extract fitted params:
                    depthRatio = Q(1);
                    impactParam(i) = Q(2);
                    transitTime = Q(3);
                    t0 = Q(4);
                    linearLimbDarkening(i) = Q(5);
                    quadLimbDarkening(i) = Q(6);
                catch e
                    % Fitting failed, go with defaults.
                    %warning(" > Fitting Failed!");
                    %warning(e.message);
                    Q = Q0;
                end
            end
            
            transitDuration(i) = transitTime;
            transitDepth(i) = 1e6 * depthRatio;
            transitCurveFit(i) = goodness_of_fit;
            % Compute Planet Parameters:
            planetRadius(i) = in.data.StarRadius(find(in.data.TIC_ID == TIC,1)) * sqrt(depthRatio) * Rsol/REarth;
        end
        
        % Mark Validity:
        validObjects(i) = valid;
    end

    % Trim data from invalid objects:
    N_peaks = N_peaks(validObjects);
    PPPeriod = PPPeriod(validObjects);
    transitDepth = transitDepth(validObjects);
    transitDuration = transitDuration(validObjects);
    transitEdgeTime = transitEdgeTime(validObjects);
    planetRadius = planetRadius(validObjects);
    impactParam = impactParam(validObjects);
    linearLimbDarkening = linearLimbDarkening(validObjects);
    quadLimbDarkening = quadLimbDarkening(validObjects);
    transitCurveFit = transitCurveFit(validObjects);

    %% Create New Dataset:
    % Columns of Output Table (name of fields following structs,
    % 'class' must be first and is not present in structs):
    columns = {'class', 'TICID', 'mag', 'Teff', 'Rstar', 'Lstar', 'Npeaks', 'PPPeriod','Depth','t_trans','t_edge', 'Rp', 'b', 'gam1','gam2', 'AR2'};
    columnDesc = {'class', 'TIC_ID', 'TESSMagnitude', 'StarTemp', 'StarRadius', 'StarLuminosity', 'NFluxRegions', 'PPPeriod','TransitDepth','TransitDuration','TransitEdgeTime', 'PlanetRadius', 'ImpactParam', 'Limb1','Limb2', 'CurveFit'}; % human-readable headers for output table
    %columns = {'class', 'TICID', 'mag', 'Teff', 'Rstar', 'Lstar', 'fit', 'depth', 'period', 'number', 'duration', 'riseTime'};
    %columnDesc = {'class', 'TIC_ID', 'TESSMagnitude', 'StarTemp', 'StarRadius', 'StarLuminosity', 'FitRSquared', 'DepthPPM', 'Period', 'Number', 'Duration', 'RiseTime'}; % human-readable headers for output table

    out.data = table;
    out.data{:,:} = in.data{validObjects,:}; % copy over valid entries
    out.data{:,end+1:end+10} = [N_peaks, PPPeriod, transitDepth, transitDuration, transitEdgeTime, planetRadius, impactParam, linearLimbDarkening,quadLimbDarkening, transitCurveFit];

    % Export Dataset:
    out.data.Properties.VariableNames = columnDesc; % Label Table
    writetable(out.data, './data/data_with_transits.csv');
end

% Fit transit curve:
% Error function to fit the curve:
function e = fit_error(qq, binned_time,binned_flux, period)
    % Compute current curve:
    fit_flux = transit_curve(binned_time.*period, qq);
    % Compute Error:
    e = sum((binned_flux - fit_flux).^2);     
end

% Fit transit curve:
% Error function to fit the curve.
% For reduced-size Q (excludes t0 and t_transit).
function e = fit_error_reduced(qq, binned_time,binned_flux, period, Q0)
    % Generate full-size (non-reduced) Q:
    qq = [qq(1); qq(2); qq(3); Q0(4); qq(4);qq(5)];
    % Compute current curve:
    fit_flux = transit_curve(binned_time.*period, qq);
    % Compute Error:
    e = sum((binned_flux - fit_flux).^2);     
end

% Error function to evaluate the curve fit (adjusted R-squared / coefficient of determination):
function AR2 = AR2(qq, binned_time,binned_flux, period)
    % Compute current curve:
    fit_flux = transit_curve(binned_time.*period, qq);
    % Compute Error:
    SS_res=sum( (binned_flux - fit_flux).^2 );
    SS_tot=sum( (binned_flux - mean(binned_flux)).^2 );

    R2 = 1 - SS_res/SS_tot;
    n = length(binned_flux);
    k = numel(qq);
    AR2 = 1 - (1 - R2) * (n-1)/(n-k-1);
end

% Constraint function on Q-space ( c(x) <= 0, ceq = 0 ).
function [c,ceq] = Q_con(qq)
    c = qq(2) ./ (1 + sqrt(qq(1))) - 1;
    ceq = 0; % ignore this constraint
end

function F = transit_curve(t, Q)
    persistent eps N_int;
    if isempty(eps)
        eps = 1e-6; % Differentiation epsilon
        N_int = 1000; % Number of integration trapezoids
    end
    
    depthRatio = Q(1);
    impactParam = Q(2);
    transitTime = Q(3); % in days (or just same units as t,t0)
    t0 = Q(4);
    gam1 = Q(5);
    gam2 = Q(6);
    
    % Convert time parameterization into z-parameterization:
    z = sqrt(4*((t-t0)./transitTime).^2 + impactParam.^2); % Both z(t) work. This is much faster (but only fits one transit dip).
    %z = (period/pi/transitTime) * sqrt( (sin(2*pi*(t-t0)/period))^2 + (pi*impactParam*transitTime/period * cos(2*pi*(t-t0)/period))^2 );
    
    p = sqrt(depthRatio);

    % Derivative of Fe(p/r,z/r)*r^2 with respect to r.
    function d = dFer2_dr(r)
        der = @(r) Fe(p./r, z./r) .* r.^2;
        if r < eps
            d = 0;
        else
            d = (der(r+eps/2) - der(r-eps/2)) / eps;
        end
    end

    % Limb-Darkening Intensity Curve:
    function I = I(r)
        I = 1 - gam1*(1-sqrt(1-r.^2)) - gam2*(1-sqrt(1-r.^2)).^2;
    end

    % Integrate simultaneously at same resolution (in a manner compatible
    % with vector input for t).
    rs = linspace(eps,1, N_int);
    int1 = zeros(numel(rs),1); % First Integrand
    int2  = zeros(numel(rs), numel(z)); % Second Integrand
    for i = 1:numel(rs)
        r = rs(i);
        int1(i) = 2 * r .* I(r);
        int2(i,:) = I(r) .* dFer2_dr(r);
    end
    
    denom = trapz(rs,int1);
    F = zeros(size(z));
    for i = 1:numel(z)
        F(i) = trapz(rs,int2(:,i)) ./ denom;
    end
end

% Undarkened / Uniform source flux curve (helper function):
function Fe = Fe(p,z)
    k0 = acos((z.^2+p.^2-1)./2./p./z);
    k1 = acos((z.^2-p.^2+1)./2./z);

    % Vector-safe piece-wise:
    lame = zeros(size(z));
    
    region = z <= (p-1);
    lame(region) = 1;
    
    region = z <= (1-p);
    if numel(p) > 1
        pregion = region;
    else
        pregion = 1;
    end
    lame(region) = p(pregion).^2;
    
    region = abs(1-p)< z & z <= (1+p);
    if numel(p) > 1
        pregion = region;
    else
        pregion = 1;
    end
    
    lame(region) = ( k0(region).*p(pregion).^2 + k1(region) - sqrt((4*z(region).^2 - (1+z(region).^2-p(pregion).^2).^2)./4 ) )./ pi;    
    
    Fe = 1 - lame;
end