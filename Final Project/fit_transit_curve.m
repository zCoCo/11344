% Used to explore the fit of a transit curve in an instance of the dataset.
% (for error analysis exploration).
function fit_transit_curve()
    % Settings:
    TIC = 149603524; % TICID of Star to Lookup

    disp("Fetching Transit Curve Data . . . ");
    
    % Load Star Data:
    starData = readtable('./data/star_data.csv');
    % Load Transit Period Table:
    transitPeriods = readtable(char("./data/transit_periods.csv"));
    
    % Fetch Relevant Metadata:
    valid = true; % Whether valid data exists for all relevant features of this object
    period = transitPeriods(transitPeriods.TICID == TIC, :).Period;
    if numel(period) > 0
        PPPeriod = string(period{1});
        PPPeriod = sscanf(PPPeriod, '%f');
    else
        PPPeriod = 0;
        valid = false; % Period not in table (periodogram couldn't be built)
    end
    
    % Load Transit Curve:
    if valid
        try
            transit = readtable(char("./transits/tr_"+TIC+".csv"));
        catch e
            valid = false;
            disp("No Valid Folded Transit Curve for TIC="+TIC);
            warning(e.message);
        end
    end
    
    if valid
        star = find(starData.TIC_ID == TIC);
        if numel(star) > 0
            Rstar = starData.StarRadius(star(1));
            Lstar = starData.StarLuminosity(star(1));
        else
            valid = false; % Period not in table (periodogram couldn't be built)
        end
    end
    
    % Compute Relevant Properties:
    if valid
        binned_flux = smooth(transit.flux, floor(numel(transit.flux)/25), 'moving');
        binned_time = transit.time;
            
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
         tr_depth = 1e6 * depthRatio;
         transitDepth = tr_depth;
        
         % Compute Timing:
        if tr_depth == 0 % No significant transit observed
            transitEdgeTime = 0;
            transitDuration = 0;
        else
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
            transitEdgeTime = 2*mean([t_fall, t_rise]);
            transitDuration = (t_50rise - t_50drop) + transitEdgeTime;
            [~, min_idx] = min(binned_flux);
            t0 = binned_time(min_idx);
        end
         
        % Stellar Mass (in SI units):
        if (Lstar/0.23)^(1/2.3) < 0.43
            Mstar = (Lstar/0.23)^(1/2.3);
        elseif 0.43 <= Lstar^0.25 && Lstar^0.25 < 2
            Mstar = Lstar^0.25;
        elseif Lstar / 32000 > 55
            Mstar = Lstar / 32000;
        else
            Mstar = (Lstar/1.4)^(1/3.5);
        end
    end
    
    % Get initial transit curve parameters:
    period = PPPeriod;
    impactParam = 0.01;
    transitTime = transitDuration*period;
    t0 = t0*period;
    Q0 = [depthRatio, impactParam, transitTime, t0, 0.5,0.1]';
    
    % Fit transit curve:
    % Error function to fit the curve:
    function e = fit_error(qq)
        % Compute current curve:
        fit_flux = transit_curve(binned_time.*period, qq);
        % Compute Error:
        e = sum((binned_flux - fit_flux).^2);     
    end

    % Fit transit curve:
    % Error function to fit the curve.
    % For reduced-size Q (excludes t0 and t_transit).
    function e = fit_error_reduced(qq)
        % Generate full-size (non-reduced) Q:
        qq = [qq(1); qq(2); qq(3); t0; qq(4);qq(5)];
        % Compute current curve:
        fit_flux = transit_curve(binned_time.*period, qq);
        % Compute Error:
        e = sum((binned_flux - fit_flux).^2);     
    end
    
    % Error function to evaluate the curve fit (adjusted R-squared / coefficient of determination):
    function AR2 = AR2(qq)
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

    disp("Fitting Transit Curve . . . ");
    
    % Get baseline goodness-of-fit:
    goodness_of_fit = AR2(Q0);
    disp(" > Pre-fit Adjusted R2: ");
    disp(goodness_of_fit);

    if(goodness_of_fit < 0.07) % <- "magic" number from tuning it / trial-and-error
        % Don't even try to optimize if initial fit is absolutely horrible (to save time).
        warning(" > Not attempting curve-fit optimization. AR2 initial is " + goodness_of_fit);
        Q = Q0;
    else
        try
            Q0_reduced = Q0([1 2 3 5 6]);
            options = optimoptions('fmincon', 'Display','iter', 'Algorithm','active-set', 'MaxIterations', 30, 'TolFun', 1e-5);
            Q = fmincon(...
                    @fit_error_reduced, ... % minimize least-squared error
                    Q0_reduced, ... % initial condition (from analyzing folded curve)
                    [0,0,0,1,1; ... % gam1 + gam2 < 1 -> Q5 + Q6 < 1
                     ...%0,0,1/period,2/period,0,0; ... % (t0-ttrans/2)/period >= -0.5 (not in reduced problem)
                    ],...%0,0,1/period,-2/period,0,0], ... % (t0+ttrans/2)/period <= 0.5 (not in reduced problem)
                    [0.999],...%; 1; 1], ...
                    [],[], ... % <- ignore lin. eq. constraints,
                    [0.75*depthRatio,  0, 0.75*transitTime, 0,0], ... % lower bound (treat as tuning, keep somewhat near recovered params)
                    [1.25*depthRatio, Inf, 1.25*transitTime, 1,1], ... % upper bound (treat as tuning, keep somewhat near recovered params).
                    @Q_con, ...                                                       % upper bound on impactParam is handled here
                    options ...
                );
            % Recreate full-Q: 
            Q = [Q(1); Q(2); Q(3); t0; Q(4); Q(5)];
            goodness_of_fit = AR2(Q);
        catch e
            % Fitting failed, go with defaults.
            warning(" > Fitting Failed!");
            warning(e.message);
            Q = Q0;
        end
    end
    disp(" > Post-fit Adjusted R2: ");
    disp(goodness_of_fit);
    
    % Extract fitted params:
    
    depthRatio = Q(1);
    impactParam = Q(2);
    transitTime = Q(3);
    t0 = Q(4);
    gam1 = Q(5);
    gam2 = Q(6);
    
    tr = @(t) transit_curve(t, Q);
    z = @(t) sqrt(4*((t-t0)/transitTime)^2 + impactParam^2);
    tr_p = @(p) tr(p*period); % parameterized by phase (period fraction)
    z_p = @(p) z(p*period); % parameterized by phase
    
    figure();
    hold on
        scatter(binned_time, binned_flux);
        ps = linspace(-0.5,0.5, 1000);
        trs = tr_p(ps);
        plot(ps,trs);
        %fplot(tr_p, [-0.5, 0.5]);
        
        t_beg = t_50drop - (t_100drop-t_50drop);
        plot([t_beg,t_beg],ylim,'k');
        plot([t_50drop,t_50drop],ylim,'k');
        plot([t_100drop,t_100drop],ylim,'k');
        plot([t_100rise,t_100rise],ylim,'k');
        plot([t_50rise,t_50rise],ylim,'k');
        t_fin = t_50rise + (t_50rise-t_100rise);
        plot([t_fin,t_fin],ylim,'k');
        
        plot(xlim,[tr_base,tr_base],'k');
        thresh1 = tr_base + 0.5*(tr_min-tr_base);
        plot(xlim,[thresh1,thresh1],'k');
        thresh2 = tr_base + 1.0*(tr_min-tr_base);
        plot(xlim,[thresh2,thresh2],'k');
        plot(xlim,[tr_min,tr_min],'k');
        
        xlabel('Phase [x$T_{\mathrm{period}}$]', 'Interpreter','latex');
        ylabel('Normalized Flux', 'Interpreter','latex');
        title(char("TIC " + TIC), 'Interpreter','latex');
        
        yyaxis right
        fplot(z_p, [-0.5, 0.5]);
        ylabel('$z=\frac{d}{R_\star}$ Phase Parameterization', 'Interpreter', 'latex');
end

function F = transit_curve(t, Q)
    persistent eps N_int;
    if isempty(eps)
        eps = 1e-6; % Differentiation epsilon
        N_int = 1000; % Number of integration trapezoids
    end
% Constants:
%     persistent Rsol Msol G;
%     Rsol = 695510e3; % [m] Radius of the Sun
%     Msol = 1.9891e30; % [kg] Mass of the Sun
%     G = 6.674e-11; % [SI units] Newton's Gravitonation Constant
%     
%     depthRatio = Q(1);
%     period = Q(2)*86400; % days -> s
%     Mstar = Q(3) * Msol; % -> SI base units
%     Rstar = Q(4) * Rsol; % -> SI base units
%     
%     a = (G * Mstar * period^2 / 4 / (pi^2))^(1/3); % Semi-Major Axis (from Kepler's Third Law)
    depthRatio = Q(1);
    impactParam = Q(2);
    transitTime = Q(3); % in days (or just same units as t,t0)
    t0 = Q(4);
    gam1 = Q(5);
    gam2 = Q(6);
%     Rstar = Rstar * Rsol; % -> SI base units
    
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

% figure();
% hold on;
% zs =  linspace(0,1, 250);
% colors = parula(numel(zs));
% for i = 1:numel(zs)
%     z = zs(i);
%     disp(z);
%     fplot(@(r) dFer2_dr(r,z), [1e-5,1], 'Color', colors(i, :));
% end
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'UpdateFcn',@myupdatefcn)
% 
% function txt = myupdatefcn(~,event_obj)
% % Customizes text of data tips
% 
% txt = char("z" + event_obj.Target.XRange(1) * 1e5);
%       
% end

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
    
%     region = abs(1-p)< z & z <= (1+p);
%     for zz = z(region)
%         
%     end
    
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

    if min(lame) < 0
        disp("break");
    end     
    
    Fe = 1 - lame;
end