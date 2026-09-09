clear; clc; close all;

%% ===================== USER SETTINGS ============================

% --- (1) Choose how to define the model: via u(n) or via f(n) directly.
%     Set MODEL_MODE to 'u' or 'f'.
MODEL_MODE = 'u';

% (1a) If MODEL_MODE = 'u': define u(n) as an anonymous function of n
%      (n = 1, 2, 3, ...). u(n) must be positive for all n you evaluate.

%      Example: the Evans "explicitly solvable" case, u(n) = beta*(1+b/n)
beta_param = 1;
b_param    = 5;
[a, b, c] = deal(3,10,44); % We need c>a+b+1
% u_fn = @(n) beta_param * (1 + b_param ./ n);
u_fn = @(n) beta_param*(n*(n+c-1))/((n+a-1)*(n+b-1));



SAVE_FIGURE = false;  % 
SHOW_TITLE = false;
EXACT_RHOC_FORMULA = 'custom';   % 'evans', 'custom', or 'none'
xmax=50;
[ymin, ymax] = deal(1e-30, 1.3);
rho = 2;                 % target density rho = N/L
L_values = [10, 30, 100, 300, 1000];  % list of L values for the canonical curves

f_fn = @(n) exp(-n .* log(beta_param) + gammaln(n + 1) - gammaln(n + b_param + 1) + gammaln(b_param + 1));


% --- (3) Numerical settings.
Nmax_pad = 1.3;   % Nmax = ceil(max(N) * Nmax_pad), extra headroom for GCE tail
Nmax_cap = 9000;  % hard cap to keep runtime/memory bounded

% --- (4) Figure appearance (fully yours to edit).
PLOT_MODE = 'scaled';  % 'linear' (p(n) vs n, log y-axis) or 'loglog' (log p(n)
                        % vs log n, both axes log-scaled -- shows power-law
                        % or 'scaled' (log p(n) vs n/N)

% lineColors = [0.20 0.55 0.85;   % smallest L
%     0.10 0.65 0.45;
%     0.85 0.70 0.15;
%     0.85 0.25 0.10];  % largest finite L


lineColors = [0.90 0.82 0.95;   % lightest purple
    0.72 0.58 0.85;
    0.55 0.35 0.70;
    0.38 0.18 0.50;
    0.20 0.05 0.28];  % near black-purple

% black = limiting/reference curve
gceColor    = [0 0 0];
gceLineWidth = 2.5;
ceLineWidth  = 1.5;
markerSize   = 5;
fontSizeAxis  = 12;
fontSizeTitle = 13;
fontSizeLegend = 10;

%% ===================== BUILD f(n) ================================

Ns_needed = round(rho .* L_values);
maxN_needed = max(Ns_needed);
Nmax = min(Nmax_cap, max(50, ceil(maxN_needed * Nmax_pad)));
if maxN_needed > Nmax_cap
    error('zrp_ce_gce:tooLarge', ...
        'N = rho*L = %d exceeds the cap of %d. Increase Nmax_cap or reduce L/rho.', ...
        maxN_needed, Nmax_cap);
end

switch MODEL_MODE
    case 'u'
        f = build_f_from_u(u_fn, Nmax);
    case 'f'
        nvec = (0:Nmax)';
        f = f_fn(nvec);
        f = f(:)'; 
        if abs(f(1) - 1) > 1e-9
            warning('zrp_ce_gce:f0', 'f(0) = %.6g, expected 1. Continuing anyway.', f(1));
        end
    otherwise
        error('zrp_ce_gce:badMode', 'MODEL_MODE must be ''u'' or ''f''.');
end

if any(~isfinite(f)) || any(f < 0)
    error('zrp_ce_gce:badF', 'f(n) produced a non-finite or negative value; check u(n)/f(n).');
end

%% ===================== CANONICAL ENSEMBLE ========================

sortedL = sort(L_values);
nSeries = numel(sortedL);
ceData = struct('L', {}, 'N', {}, 'n', {}, 'p', {}, 'meanUCE', {});

fprintf('Computing canonical marginals via repeated squaring...\n');
for i = 1:nSeries
    L = sortedL(i);
    N = round(rho * L);
    tic;
    [p, meanUCE] = pCE_repeated_squaring(f, L, N);
    dt = toc;
    fprintf('  L = %6d, N = %5d : %.3f s\n', L, N, dt);
    ceData(i).L = L;
    ceData(i).N = N;
    ceData(i).n = (0:N)';
    ceData(i).p = p(:);
    ceData(i).meanUCE = meanUCE;
end

%% ===================== GRAND-CANONICAL ENSEMBLE ==================

[z_star, beta_est, rho_at_beta] = solve_fugacity(f, rho);

switch EXACT_RHOC_FORMULA
    case 'evans'
        if b_param > 2
            rho_c_display = 1 / (b_param - 2);
        else
            rho_c_display = Inf;
        end
        rho_c_isExact = true;
    case 'custom'
        % u(n) = beta*n*(n+c-1) / ((n+a-1)*(n+b-1)); condensation iff c > a+b+1
        if c > a + b + 1
            rho_c_display = a*b / (c - a - b - 1);
        else
            rho_c_display = Inf;
        end
        rho_c_isExact = true;
    otherwise
        rho_c_display = rho_at_beta;
        rho_c_isExact = false;
end

hasGCE = ~isnan(z_star);
if hasGCE
    plotMaxN = max(Ns_needed);
    extentN = min(Nmax, ceil(plotMaxN * 1.25) + 2);
    n_gce = (0:extentN)';
    p_gce = pGCE_from_z(f, z_star, extentN);
    fprintf('GCE fugacity z(rho=%.4g) = %.6f  (beta ~ %.4f)\n', rho, z_star, beta_est);
else
    fprintf(['rho = %.4g exceeds rho_c (approx density at beta = %.4g): ' ...
             'no fluid-phase GCE solution -- system is in the condensed regime.\n'], ...
             rho, rho_at_beta);
end

%% ===================== FIGURE =====================================

figure('Color', 'w', 'Position', [100 100 800 600]);
hold on; box on;

legendEntries = {};
legendHandles = [];

switch PLOT_MODE
    case 'linear'
        
        if hasGCE
            h_gce = plot(n_gce, p_gce, '-', 'Color', gceColor, 'LineWidth', gceLineWidth);
            legendHandles(end+1) = h_gce; 
            legendEntries{end+1} = 'GCE: p(n; z^*)';
        end

        for i = 1:nSeries
            col = lineColors(mod(i-1, size(lineColors,1)) + 1, :);
            nz = ceData(i).p > 0;
            xline_ = ceData(i).n(nz);
            yline_ = ceData(i).p(nz);
            h = plot(xline_, yline_, '-o', 'Color', col, 'LineWidth', ceLineWidth, ...
                'MarkerFaceColor', col, 'MarkerSize', markerSize);

            if ~isempty(yline_)
                plot([ceData(i).N ceData(i).N], [yline_(end), yline_(end)*1e-3], '--', ...
                    'Color', col, 'LineWidth', ceLineWidth, 'HandleVisibility', 'off');
            end
            legendHandles(end+1) = h;
            legendEntries{end+1} = sprintf('CE: L = %d, N = %d', ceData(i).L, ceData(i).N);
        end

        set(gca, 'YScale', 'log', 'XScale', 'linear');
        xlabel('n', 'FontSize', fontSizeAxis);
        ylabel('p(n)', 'FontSize', fontSizeAxis);

        
        allP = [];
        if hasGCE, allP = [allP; p_gce(p_gce > 0)]; end
        for i = 1:nSeries
            allP = [allP; ceData(i).p(ceData(i).p > 0)]; 
        end
        if ~isempty(allP)
            yFloor = min(allP);
            ylim([10^(floor(log10(yFloor)) - 0.3), 1.3]);
        end
    %%%%%%%%%
    case 'loglog'
        
        if hasGCE
            nz = (n_gce >= 1) & (p_gce > 0);
            h_gce = plot(n_gce(nz), p_gce(nz), '-', 'Color', gceColor, 'LineWidth', gceLineWidth);
            legendHandles(end+1) = h_gce;
            legendEntries{end+1} = 'GCE: p(n; z^*)';
        end

        for i = 1:nSeries
            col = lineColors(mod(i-1, size(lineColors,1)) + 1, :);
            nz = (ceData(i).n >= 1) & (ceData(i).p > 0);
            xline_ = ceData(i).n(nz);
            yline_ = ceData(i).p(nz);
            h = plot(xline_, yline_, '-o', 'Color', col, 'LineWidth', ceLineWidth, ...
                'MarkerFaceColor', col, 'MarkerSize', markerSize);
            legendHandles(end+1) = h;
            legendEntries{end+1} = sprintf('CE: L = %d, N = %d', ceData(i).L, ceData(i).N);
        end

        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('n', 'FontSize', fontSizeAxis);
        ylabel('p(n)', 'FontSize', fontSizeAxis);
    %%%%%%% 
    case 'scaled'
        
        for i = 1:nSeries
            col = lineColors(mod(i-1, size(lineColors,1)) + 1, :);
            nz = ceData(i).p > 0;
            xline_ = ceData(i).n(nz) / ceData(i).N;
            yline_ = ceData(i).p(nz);
            h = plot(xline_, yline_, '-o', 'Color', col, 'LineWidth', ceLineWidth, ...
                'MarkerFaceColor', col, 'MarkerSize', markerSize);
            legendHandles(end+1) = h;
            legendEntries{end+1} = sprintf('CE: L = %d, N = %d', ceData(i).L, ceData(i).N);
        end
        
        
        if rho_c_isExact && isfinite(rho_c_display) && rho > rho_c_display
            peakLoc = 1 - rho_c_display / rho;
            xline(peakLoc, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
        
        set(gca, 'YScale', 'log', 'XScale', 'linear');
        xlabel('n / N', 'FontSize', fontSizeAxis);
        ylabel('p(n)', 'FontSize', fontSizeAxis);

        allP = [];
        for i = 1:nSeries
            allP = [allP; ceData(i).p(ceData(i).p > 0)]; 
        end
        if ~isempty(allP)
            yFloor = min(allP);
            ylim([10^(floor(log10(yFloor)) - 0.3), 1.3]);
        end
        
    otherwise
        error('zrp_ce_gce:badPlotMode', 'PLOT_MODE must be ''linear'', ''loglog'', or ''scaled''.');
end


if rho_c_isExact
    rhoCTexSymbol = '=';
else
    rhoCTexSymbol = '\approx';
end
if rho_c_isExact && isinf(rho_c_display)
    rhoCValueStr = '\infty';
else
    rhoCValueStr = sprintf('%.4g', rho_c_display);
end

if hasGCE
    titleStr = sprintf('ZRP single-site marginal: u(n) = 1+%.4g/n, \\rho = %.4g,  \\rho_c %s %s', ...
    b_param, rho, rhoCTexSymbol, rhoCValueStr);
else
    titleStr = sprintf('ZRP single-site marginal (CONDENSED): u(n)=1+%.4g/n, \\rho = %.4g > \\rho_c %s %s', ...
    b_param, rho, rhoCTexSymbol, rhoCValueStr);
end

if SHOW_TITLE
    title(titleStr, 'FontSize', fontSizeTitle, 'Interpreter', 'tex');
end
legend(legendHandles, legendEntries, 'Location', 'northeast', ...
'FontSize', fontSizeLegend, 'Box', 'off');

grid off;
set(gca, 'FontSize', fontSizeAxis, 'GridAlpha', 0.15);
hold off;

%%%% Setting x-range
if strcmp(PLOT_MODE, 'scaled')
    xlim([0 1.05]);   
else
    xlim([0 xmax]);
end
ylim([ymin ymax]); % set y-axis range 

fprintf('\nDone. Adjust lineColors / fontSize* / figure Position above for further control.\n');
fprintf('Switch PLOT_MODE to ''loglog'' to see power-law decay / condensate upturn.\n');
fprintf('Save with, e.g.:  print(gcf, ''zrp_figure.pdf'', ''-dpdf'', ''-bestfit'')\n');
if SAVE_FIGURE
    saveFolder = '/Users/dave/Documents/imperial_applied_msc/dissertation!/dissertation matlab/dissertation_figures';
    figName = sprintf('b%g_rho%g_%s', b_param, rho, PLOT_MODE);
    print(gcf, fullfile(saveFolder, [figName '.png']), '-dpng', '-r300');
end

%% Printing the caption for LaTeX
if rho_c_isExact
    if isinf(rho_c_display)
        rhoCCaptionStr = '\infty';
    else
        rhoCCaptionStr = sprintf('%.4g', rho_c_display);
    end
    rhoCCaptionRelation = '=';
else
    rhoCCaptionStr = sprintf('%.4g', rho_c_display);
    rhoCCaptionRelation = '\approx';
end

if rho < rho_c_display
    rhoRelationWord = '<';
else
    rhoRelationWord = '>';
end

uDescr = sprintf('u(n) = 1 + %.4g/n', b_param); 

LNPairs = cell(1, nSeries);
for i = 1:nSeries
    LNPairs{i} = sprintf('(%d,%d)', ceData(i).L, ceData(i).N);
end
LNPairsStr = strjoin(LNPairs, ', ');

captionStr = sprintf(['Single-site marginal $p(n)$ for the ZRP with jump rate ' ...
    '$%s$ at density $\\rho = %.4g %s \\rho_c %s %s$. ' ...
    'Black: GCE marginal $p(n; z^*)$. Coloured: CE marginals from exact ' ...
    'enumeration at $(L,N) = %s$. Scale: log-linear plot.'], ...
    uDescr, rho, rhoRelationWord, rhoCCaptionRelation, rhoCCaptionStr, LNPairsStr);

fprintf('\n--- LaTeX caption (paste into \\caption{...} in Overleaf) ---\n\n');
fprintf('%s\n\n', captionStr);
