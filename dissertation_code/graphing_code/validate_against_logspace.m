%% validate_against_logspace.m
%
% Cross-checks pCE_repeated_squaring (the fast method used throughout this
% toolkit) against an independent, deliberately slow, straightforward
% log-space direct recursion for the canonical partition function. This
% is the same style of validation used during development of the
% companion JavaScript/React tool, and is included here so you can
% independently confirm correctness on your own machine before trusting
% the fast method for your dissertation figures.
%
% Run this once after downloading the toolkit; it should print relative
% differences on the order of 1e-13 or smaller across several (L, N)
% combinations, and should NOT error.

clear; clc;

u_fn = @(n) 1 * (1 + 5 ./ n); % Evans exact case, beta=1, b=5, rho_c = 1/3

testCases = [10, 2; 30, 6; 100, 20; 300, 60];

fprintf('%-8s %-8s %-16s %-14s\n', 'L', 'N', 'max rel diff', 'sum(p_fast)');
for i = 1:size(testCases, 1)
    L = testCases(i, 1);
    N = testCases(i, 2);

    f = build_f_from_u(u_fn, N);
    [p_fast, ~] = pCE_repeated_squaring(f, L, N);

    logf = build_logf_slow(u_fn, N);
    table = compute_logZ_table_slow(logf, L, N);
    p_trusted = pCE_from_table_slow(logf, table, L, N);

    maxRelDiff = 0;
    for n = 0:N
        if p_trusted(n + 1) > 0
            rd = abs(p_trusted(n + 1) - p_fast(n + 1)) / p_trusted(n + 1);
            maxRelDiff = max(maxRelDiff, rd);
        end
    end

    fprintf('%-8d %-8d %-16.3e %-14.12f\n', L, N, maxRelDiff, sum(p_fast));

    if maxRelDiff > 1e-8
        warning('validate_against_logspace:mismatch', ...
            'Relative difference %.3e at L=%d, N=%d is larger than expected.', ...
            maxRelDiff, L, N);
    end
end

fprintf('\nIf all "max rel diff" values above are ~1e-12 or smaller, the fast\n');
fprintf('method is working correctly on your system.\n');

%% ---- Slow, independent reference implementation (log-space recursion) ----

function logf = build_logf_slow(u_fn, Nmax)
    logf = -inf(1, Nmax + 1);
    logf(1) = 0;
    acc = 0;
    for n = 1:Nmax
        acc = acc - log(u_fn(n));
        logf(n + 1) = acc;
    end
end

function table = compute_logZ_table_slow(logf, L, Nmax)
    prev = -inf(1, Nmax + 1);
    prev(1) = 0;
    table = cell(1, L + 1);
    table{1} = prev;
    for l = 1:L
        cur = -inf(1, Nmax + 1);
        for n = 0:Nmax
            terms = [];
            for k = 0:n
                if logf(k + 1) == -inf || prev(n - k + 1) == -inf
                    continue;
                end
                terms(end + 1) = logf(k + 1) + prev(n - k + 1); %#ok<AGROW>
            end
            cur(n + 1) = logsumexp_slow(terms);
        end
        table{l + 1} = cur;
        prev = cur;
    end
end

function p = pCE_from_table_slow(logf, table, L, N)
    logZL = table{L + 1}(N + 1);
    rowLm1 = table{L};
    p = zeros(N + 1, 1);
    if logZL == -inf
        return;
    end
    for n = 0:N
        idx = N - n;
        if idx < 0 || idx + 1 > numel(rowLm1) || rowLm1(idx + 1) == -inf || logf(n + 1) == -inf
            continue;
        end
        p(n + 1) = exp(logf(n + 1) + rowLm1(idx + 1) - logZL);
    end
end

function s = logsumexp_slow(x)
    if isempty(x)
        s = -inf;
        return;
    end
    m = max(x);
    if ~isfinite(m)
        s = -inf;
        return;
    end
    s = m + log(sum(exp(x - m)));
end
