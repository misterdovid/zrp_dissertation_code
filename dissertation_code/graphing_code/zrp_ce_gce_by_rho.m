%% zrp_ce_gce_by_rho.m
%
% Zero-Range Process: canonical (CE) single-site marginal p_CE(n) at a
% FIXED (large) L, shown for SEVERAL densities rho -- including, if you
% like, densities above rho_c (the condensed regime; those series are
% still plotted and are flagged as condensed, not excluded). A SINGLE
% grand-canonical reference curve is also shown: the GCE marginal
% evaluated AT CRITICALITY, p_GCE(n; beta), i.e. at the largest allowed
% fugacity z = beta (only defined when rho_c is finite).
%
% This is the natural counterpart to zrp_ce_gce.m, which instead fixes rho
% and varies L: here we fix L (large) and vary rho, showing how the CE
% marginal changes across a range of densities you choose -- fluid phase,
% near criticality, or condensed -- relative to the single critical GCE
% curve.
%
% Works for ANY b_param (including b_param <= 2, where rho_c = Inf, there
% is no condensation at any density, and no critical GCE curve is drawn)
% as well as any rho you specify (including rho >= rho_c when rho_c is
% finite).
%
% ------------------------------------------------------------------
% BACKGROUND (see Evans & Hanney, J. Phys. A 38 (2005) R195)
%
%   f(0) = 1,  f(n) = prod_{k=1}^n 1/u(k)              for n > 0
%
%   Canonical single-site marginal (fixed L, N = round(rho*L)):
%       p_CE(n) = f(n) * Z_{L-1,N-n} / Z_{L,N}
%
%   Grand-canonical marginal at fugacity z:
%       p_GCE(n; z) = z^n f(n) / F(z)
%
%   For u(n) = beta_param*(1 + b_param/n), beta = beta_param exactly
%   (since u(n) -> beta_param as n -> infinity), and
%       rho_c = 1/(b_param - 2)   if b_param > 2,   else rho_c = Inf.
%   (Evans & Hanney, Eq. 59.)
%
% ------------------------------------------------------------------
% ALGORITHM / PERFORMANCE NOTE
%
% Uses the same fast, numerically-validated repeated-squaring method as
% zrp_ce_gce.m for the CE marginals (see that file's header and
% poly_pow_trunc_scaled.m for why plain O(N^2) convolution is used instead
% of FFT). The critical GCE curve is evaluated directly at the EXACT
% beta = beta_param (not a numerical estimate, and not via bisection at
% all -- there is nothing to solve for here, since criticality means
% z = beta by definition).
%
% ------------------------------------------------------------------

clear; clc; close all;

%% ===================== USER SETTINGS ============================

% --- (1) Model: MUST be the Evans form u(n) = beta_param*(1+b_param/n)
%     for this script (needed so beta, and rho_c's finiteness, are known
%     exactly rather than estimated). b_param <= 2 is fine (rho_c = Inf,
%     no condensation at any density) -- there is no restriction on
%     b_param here.
beta_param = 1;
b_param    = 5;
u_fn = @(n) beta_param * (1 + b_param ./ n);



beta_exact = beta_param;
if b_param > 2
    rho_c_exact = 1 / (b_param - 2);
    rhoCIsFinite = true;
else
    rho_c_exact = Inf;
    rhoCIsFinite = false;
end

% --- (2) Fixed (large) L, and the list of densities to plot. Any values
%     are allowed, including rho >= rho_c (condensed regime, when rho_c
%     is finite) -- those series will show a CE curve but no GCE curve,
%     and this is reported in the console output, not treated as an error.
L = 1000;
rho_values = [0.25, 1, 4];

% --- (3) Numerical settings.
Nmax_pad = 1.3;   % extra headroom above the largest N needed, for the GCE tail
Nmax_cap = 20000; % hard cap to keep runtime/memory bounded

% --- (4) Figure appearance (fully yours to edit).
PLOT_MODE = 'loglog';   % 'linear' (p(n) vs n, log y-axis) or 'loglog'
SHOW_TITLE = false;     % false for thesis figures (info goes in \caption instead) or 'scaled'
SHOW_CRITICAL_GCE = false;   % set false to hide the critical GCE curve entirely
SAVE_FIGURE = false;    % true to write the figure to disk (see save block below)
saveFolder = '/Users/dave/Documents/imperial_applied_msc/dissertation!/dissertation matlab/dissertation_figures 2';

xmax=1e4;
[ymin, ymax] = deal(1e-13, 1.3);
% Blue -> red gradient, ordered to match increasing density (first rho =
% most blue, last rho = most red). nColors automatically matches however
% many series you're plotting.
% Light blue -> pale/neutral -> light red, ordered to match increasing
% density. Interpolates through a light grey midpoint instead of purple.
nColors = numel(rho_values);
t = linspace(0, 1, nColors)';           % 0 = first (lowest rho), 1 = last (highest)
lightBlue = [0.55 0.70 0.95];
paleMid   = [0.85 0.85 0.82];
lightRed  = [0.95 0.55 0.55];

lineColors = zeros(nColors, 3);
for k = 1:nColors
    if t(k) <= 0.5
        w = t(k) / 0.5;                 % 0..1 across first half
        lineColors(k,:) = (1-w)*lightBlue + w*paleMid;
    else
        w = (t(k) - 0.5) / 0.5;         % 0..1 across second half
        lineColors(k,:) = (1-w)*paleMid + w*lightRed;
    end
end
gceColor     = [0 0 0];
gceLineWidth = 2.5;
ceLineWidth  = 1.5;
markerSize   = 5;
fontSizeAxis   = 12;
fontSizeTitle  = 13;
fontSizeLegend = 10;

%% ===================== BUILD f(n) ================================

Ns_needed = round(rho_values .* L);
maxN_needed = max(Ns_needed);
Nmax = min(Nmax_cap, max(50, ceil(maxN_needed * Nmax_pad)));
if maxN_needed > Nmax_cap
    error('zrp_ce_gce_by_rho:tooLarge', ...
        'N = rho*L = %d exceeds the cap of %d for some rho value. Increase Nmax_cap, reduce L, or use smaller rho.', ...
        maxN_needed, Nmax_cap);
end

f = build_f_from_u(u_fn, Nmax);
if any(~isfinite(f)) || any(f < 0)
    error('zrp_ce_gce_by_rho:badF', 'f(n) produced a non-finite or negative value; check u(n).');
end

nSeries = numel(rho_values);

%% ===================== CANONICAL ENSEMBLE, PER RHO =================

data = struct('rho', {}, 'N', {}, 'n', {}, 'p', {}, 'meanUCE', {});

fprintf('Computing marginals at fixed L = %d for %d density value(s)...\n', L, nSeries);
if rhoCIsFinite
    fprintf('rho_c = %.6g (exact)\n', rho_c_exact);
else
    fprintf('rho_c = Inf (b_param = %.4g <= 2: no condensation at any density)\n', b_param);
end

for i = 1:nSeries
    rho_i = rho_values(i);
    N = round(rho_i * L);

    tic;
    [p, meanUCE] = pCE_repeated_squaring(f, L, N);
    dtCE = toc;

    if rhoCIsFinite && rho_i >= rho_c_exact
        fprintf('  rho = %6.4g, N = %6d : CE %.3fs (condensed: rho >= rho_c = %.6g)\n', ...
            rho_i, N, dtCE, rho_c_exact);
    else
        fprintf('  rho = %6.4g, N = %6d : CE %.3fs\n', rho_i, N, dtCE);
    end

    data(i).rho = rho_i;
    data(i).N = N;
    data(i).n = (0:N)';
    data(i).p = p(:);
    data(i).meanUCE = meanUCE;
end

% Single critical GCE reference curve, p_GCE(n; beta), shown once
% regardless of which/how many rho values are plotted above. Only defined
% when rho_c is finite (b_param > 2) -- with rho_c = Inf there is no
% finite critical fugacity to evaluate at.
hasCriticalGCE = rhoCIsFinite && SHOW_CRITICAL_GCE;
if hasCriticalGCE
    plotMaxN = max(Ns_needed);
    extentN = min(Nmax, ceil(plotMaxN * 1.25) + 2);
    n_gce = (0:extentN)';
    p_gce = pGCE_from_z(f, beta_exact, extentN);
    fprintf('Critical GCE: z = beta = %.6g exactly (rho_c = %.6g)\n', beta_exact, rho_c_exact);
    fprintf('sum(p_gce, truncated at N=%d) = %.10f  (should be close to 1; if noticeably\n', extentN, sum(p_gce));
    fprintf('  less than 1, increase Nmax_pad/Nmax_cap -- the power-law tail at\n');
    fprintf('  criticality decays slowly and needs more terms to sum accurately)\n');
end

%% ===================== FIGURE =====================================

figure('Color', 'w', 'Position', [100 100 800 600]);
hold on; box on;

legendEntries = {};
legendHandles = [];

switch PLOT_MODE
    case 'linear'
        if hasCriticalGCE
            h_gce = plot(n_gce, p_gce, '-', 'Color', gceColor, 'LineWidth', gceLineWidth); %#ok<*SAGROW>
            legendHandles(end+1) = h_gce;
            legendEntries{end+1} = 'GCE at criticality: p(n; \beta)';
        end

        for i = 1:nSeries
            col = lineColors(mod(i-1, size(lineColors,1)) + 1, :);
            nz = data(i).p > 0;
            xline_ = data(i).n(nz);
            yline_ = data(i).p(nz);
            h = plot(xline_, yline_, '--o', 'Color', col, 'LineWidth', ceLineWidth, ...
                'MarkerFaceColor', col, 'MarkerSize', markerSize);
            if ~isempty(yline_)
                plot([data(i).N data(i).N], [yline_(end), yline_(end)*1e-3], '--', ...
                    'Color', col, 'LineWidth', ceLineWidth, 'HandleVisibility', 'off');
            end
            legendHandles(end+1) = h;
            if rhoCIsFinite && data(i).rho >= rho_c_exact
                legendEntries{end+1} = sprintf('CE: \\rho = %.3g, N = %d (condensed)', data(i).rho, data(i).N);
            else
                legendEntries{end+1} = sprintf('CE: \\rho = %.3g, N = %d', data(i).rho, data(i).N);
            end
        end

        set(gca, 'YScale', 'log', 'XScale', 'linear');
        xlabel('n', 'FontSize', fontSizeAxis);
        ylabel('p(n)', 'FontSize', fontSizeAxis);

        allP = [];
        if hasCriticalGCE, allP = [allP; p_gce(p_gce > 0)]; end
        for i = 1:nSeries
            allP = [allP; data(i).p(data(i).p > 0)]; %#ok<AGROW>
        end
        if ~isempty(allP)
            yFloor = min(allP);
            ylim([10^(floor(log10(yFloor)) - 0.3), 1.3]);
        end

    case 'loglog'
        if hasCriticalGCE
            nzg = (n_gce >= 1) & (p_gce > 0);
            h_gce = plot(n_gce(nzg), p_gce(nzg), '-', 'Color', gceColor, 'LineWidth', gceLineWidth);
            legendHandles(end+1) = h_gce;
            legendEntries{end+1} = 'GCE at criticality: p(n; \beta)';
        end

        for i = 1:nSeries
            col = lineColors(mod(i-1, size(lineColors,1)) + 1, :);
            nz = (data(i).n >= 1) & (data(i).p > 0);
            xline_ = data(i).n(nz);
            yline_ = data(i).p(nz);
            h = plot(xline_, yline_, '-o', 'Color', col, 'LineWidth', ceLineWidth, ...
                'MarkerFaceColor', col, 'MarkerSize', markerSize);
            legendHandles(end+1) = h;
            if rhoCIsFinite && data(i).rho >= rho_c_exact
                legendEntries{end+1} = sprintf('CE: \\rho = %.3g, N = %d (condensed)', data(i).rho, data(i).N);
            else
                legendEntries{end+1} = sprintf('CE: \\rho = %.3g, N = %d', data(i).rho, data(i).N);
            end
        end

        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('n', 'FontSize', fontSizeAxis);
        ylabel('p(n)', 'FontSize', fontSizeAxis);
    case 'scaled'
        % ---- log p(n) vs n/N: rescales each CE curve onto [0,1] by its
        % own N, making the different-rho curves directly comparable. No
        % GCE curve shown here (the critical GCE has no single N to
        % rescale by that would be meaningful across all series). ----
        for i = 1:nSeries
            col = lineColors(mod(i-1, size(lineColors,1)) + 1, :);
            nz = data(i).p > 0;
            xline_ = data(i).n(nz) / data(i).N;
            yline_ = data(i).p(nz);
            h = plot(xline_, yline_, '-o', 'Color', col, 'LineWidth', ceLineWidth, ...
                'MarkerFaceColor', col, 'MarkerSize', markerSize);
            legendHandles(end+1) = h;
            if rhoCIsFinite && data(i).rho >= rho_c_exact
                legendEntries{end+1} = sprintf('CE: \\rho = %.3g, N = %d (condensed)', data(i).rho, data(i).N);
            else
                legendEntries{end+1} = sprintf('CE: \\rho = %.3g, N = %d', data(i).rho, data(i).N);
            end
        end

        set(gca, 'YScale', 'log', 'XScale', 'linear');
        xlabel('n / N', 'FontSize', fontSizeAxis);
        ylabel('p(n)', 'FontSize', fontSizeAxis);

        allP = [];
        for i = 1:nSeries
            allP = [allP; data(i).p(data(i).p > 0)]; %#ok<AGROW>
        end
        if ~isempty(allP)
            yFloor = min(allP);
            ylim([10^(floor(log10(yFloor)) - 0.3), 1.3]);
        end
    otherwise
                error('zrp_ce_gce_by_rho:badPlotMode', 'PLOT_MODE must be ''linear'', ''loglog'', or ''scaled''.');
end

if SHOW_TITLE
    if rhoCIsFinite
        titleStr = sprintf('ZRP single-site marginal at fixed L = %d: \\rho_c = %.4g', L, rho_c_exact);
    else
        titleStr = sprintf('ZRP single-site marginal at fixed L = %d: \\rho_c = \\infty', L);
    end
    title(titleStr, 'FontSize', fontSizeTitle, 'Interpreter', 'tex');
end

legend(legendHandles, legendEntries, 'Location', 'northeast', ...
    'FontSize', fontSizeLegend, 'Box', 'off');

xlim([0, xmax]);
ylim([ymin, ymax]);


grid off;
set(gca, 'FontSize', fontSizeAxis, 'GridAlpha', 0.15);
hold off;

%% ===================== LATEX CAPTION (for thesis use) ==============
%
% Builds a caption string analogous to zrp_ce_gce.m's, for this
% fixed-L / varying-rho comparison. Paste the printed LaTeX directly into
% a \caption{...} in Overleaf.

if beta_param == 1
    uDescr = sprintf('u(n) = 1 + %.4g/n', b_param);
else
    uDescr = sprintf('u(n) = %.4g(1 + %.4g/n)', beta_param, b_param);
end

if rhoCIsFinite
    rhoCCaptionStr = sprintf('\\rho_c = %.4g', rho_c_exact);
else
    rhoCCaptionStr = '\rho_c = \infty';
end

rhoNPairs = cell(1, nSeries);
anyCondensed = false;
for i = 1:nSeries
    isCondensed = rhoCIsFinite && data(i).rho >= rho_c_exact;
    if isCondensed
        rhoNPairs{i} = sprintf('(%.3g,%d)^*', data(i).rho, data(i).N); % flag condensed with *
        anyCondensed = true;
    else
        rhoNPairs{i} = sprintf('(%.3g,%d)', data(i).rho, data(i).N);
    end
end
rhoNPairsStr = strjoin(rhoNPairs, ', ');

if anyCondensed
    condensedNote = ' ($^*$ denotes a condensed, i.e. $\rho \geq \rho_c$, density)';
else
    condensedNote = '';
end

if hasCriticalGCE
    gceCaptionSentence = 'Black: GCE marginal at criticality, $p(n; \beta)$. ';
else
    gceCaptionSentence = '';
end

captionStr = sprintf(['Single-site marginal $p(n)$ for the ZRP with jump rate ' ...
    '$%s$ at fixed $L = %d$ ($%s$), for several densities. ' ...
    '%s' ...
    'Coloured: CE marginals from exact enumeration at $(\\rho, N) = %s$%s.'], ...
    uDescr, L, rhoCCaptionStr, gceCaptionSentence, rhoNPairsStr, condensedNote);

fprintf('\n--- LaTeX caption (paste into \\caption{...} in Overleaf) ---\n\n');
fprintf('%s\n\n', captionStr);
fprintf('--------------------------------------------------------------\n');

%% ===================== SAVE FIGURE (optional) =======================

if SAVE_FIGURE
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end
    rhoTag = strjoin(arrayfun(@(r) sprintf('%.3g', r), rho_values, 'UniformOutput', false), '-');
    figName = sprintf('L%d_rho%s_%s', L, rhoTag, PLOT_MODE);
    outPath = fullfile(saveFolder, [figName '.png']);
    print(gcf, outPath, '-dpng', '-r300');
    fprintf('Saved figure to %s\n', outPath);
else
    fprintf('SAVE_FIGURE is false -- figure not saved (set to true to save this run).\n');
end

fprintf('\nDone. Adjust lineColors / fontSize* / figure Position above for further control.\n');
fprintf('Switch PLOT_MODE to ''loglog'' to see power-law decay / condensate upturn more clearly.\n');
