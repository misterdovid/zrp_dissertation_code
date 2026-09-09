function [z, betaEst, rhoAtBeta] = solve_fugacity(f, rho, exactBeta, rhoCIsFinite)
%SOLVE_FUGACITY  Solve rho = z F'(z) / F(z) for the fugacity z in (0, beta),
%   where F(z) = sum f(n) z^n and beta is F's radius of convergence,
%   via bisection (rho(z) is strictly increasing on (0, beta), a standard
%   result for this problem -- see the accompanying notes).
%
%   f            : vector of single-site weights f(0), f(1), ..., f(Nmax)
%   rho          : target density
%   exactBeta    : OPTIONAL. If you know beta exactly (e.g. the Evans case
%                  u(n) = beta_param*(1+b/n), where beta = beta_param
%                  exactly), pass it here instead of relying on the
%                  numerical Aitken-extrapolated estimate. Omit (or pass
%                  []) to fall back to the numerical estimate, the only
%                  option for a general u(n)/f(n).
%   rhoCIsFinite : OPTIONAL (default true if exactBeta given, else
%                  irrelevant). Set to FALSE whenever you know rho_c is
%                  actually infinite (e.g. Evans with b<=2), even though
%                  beta itself is finite. This matters because: the
%                  density at z=beta is computed from a TRUNCATED sum (up
%                  to whatever Nmax the f(n) vector covers), so when the
%                  true rho_c is infinite, that truncated sum still
%                  produces some finite-looking number (e.g. ~60) purely
%                  as an artifact of the truncation -- NOT a real
%                  supercritical boundary. If rhoCIsFinite is false, this
%                  density-at-beta value is never used as a hard cutoff:
%                  bisection is allowed to push z arbitrarily close to
%                  beta to reach any requested rho, however large.
%
%   Returns:
%   z          : the fugacity solving the density equation, or NaN if
%                rho >= rho_c (only possible when rho_c is finite, i.e.
%                rhoCIsFinite was true or omitted with a finite exactBeta)
%   betaEst    : beta (exact, if exactBeta was supplied; else estimated)
%   rhoAtBeta  : the density attained at (or very near) z=beta -- if
%                rhoCIsFinite is false, this is a TRUNCATED-SUM ESTIMATE
%                only, not the true (infinite) rho_c, and should not be
%                reported to the user as "the" critical density.
%
%   Beta is estimated from f(n) directly via the ratio test (f(n-1)/f(n)
%   -> beta as n -> infinity for well-behaved f), using Aitken
%   extrapolation on the ratio sequence at large n to accelerate
%   convergence and correct for slow (O(1/n)) approach to the limit,
%   which a naive single-ratio estimate underestimates badly.

    haveExactBeta = nargin >= 3 && ~isempty(exactBeta) && isfinite(exactBeta);
    if nargin < 4 || isempty(rhoCIsFinite)
        rhoCIsFinite = haveExactBeta; % sensible default: if you gave an exact
                                       % beta but didn't say otherwise, assume
                                       % you mean a genuine finite rho_c
    end

    Nmax = numel(f) - 1;
    if haveExactBeta
        betaEst = exactBeta;
    else
        betaEst = estimate_beta(f);
    end

    if haveExactBeta && rhoCIsFinite
        % Check directly at z = beta first: avoids the systematic
        % undershoot from bisecting only up to beta*(1-1e-6), which
        % matters specifically for densities close to (or, for a
        % supercritical sweep, at/above) the true finite rho_c.
        dAtExactBeta = density_at_z(f, exactBeta, Nmax);
        if isfinite(dAtExactBeta)
            rhoAtBeta = dAtExactBeta;
            if rho >= dAtExactBeta
                z = NaN; % at/above the true critical density: condensed regime
                return;
            end
            % else rho < rho_c: fall through to ordinary bisection, fine
            % away from the boundary.
        end
    end

    hiCap = betaEst * (1 - 1e-6);
    if ~isfinite(hiCap) || hiCap <= 0
        hiCap = 1e6; % fallback if beta could not be estimated (e.g. beta = Inf)
    end

    lo = 1e-10;
    hi = hiCap;
    rhoAtBeta = density_at_z(f, hi, Nmax);

    if ~rhoCIsFinite
        % rho_c is known to be infinite: the truncated-sum "density at
        % hi" is NOT a real physical bound, so never treat it as one.
        % Instead, push hi progressively closer to beta until the
        % (truncated) density at hi exceeds the target rho, or we run out
        % of safe headroom -- this lets bisection reach large rho values
        % correctly rather than being capped by a truncation artifact.
        attempts = 0;
        while (~isfinite(rhoAtBeta) || rhoAtBeta < rho) && attempts < 60
            hi = hi + (betaEst - hi) * 0.5; % halve the remaining gap to beta
            if (betaEst - hi) < eps(betaEst) * 10
                break; % as close to beta as double precision allows
            end
            rhoAtBeta = density_at_z(f, hi, Nmax);
            attempts = attempts + 1;
        end
        if ~isfinite(rhoAtBeta) || rhoAtBeta < rho
            warning('solve_fugacity:truncationLimited', ...
                ['Could not reach the requested density rho = %.4g even at ' ...
                 'z very close to beta -- this is a truncated-sum limitation ' ...
                 '(increase Nmax) rather than a true condensation transition, ' ...
                 'since rho_c is infinite here. Returning the closest ' ...
                 'achievable fugacity instead.'], rho);
            z = hi;
            return;
        end
    else
        if ~isfinite(rhoAtBeta) || rhoAtBeta < rho
            z = NaN;
            return;
        end
    end

    rhoAtLo = density_at_z(f, lo, Nmax);
    if rhoAtLo >= rho
        z = lo;
        return;
    end

    for iter = 1:200
        mid = 0.5 * (lo + hi);
        d = density_at_z(f, mid, Nmax);
        if ~isfinite(d)
            hi = mid;
            continue;
        end
        if d < rho
            lo = mid;
        else
            hi = mid;
        end
    end
    z = 0.5 * (lo + hi);
end

function beta = estimate_beta(f)
%ESTIMATE_BETA  Estimate the radius of convergence of F(z) = sum f(n) z^n
%   via the ratio test, beta = lim f(n-1)/f(n), accelerated with Aitken
%   extrapolation at two different large-n scales.
    Nmax = numel(f) - 1;
    n0 = min(200, max(10, floor(Nmax / 4)));
    a1 = try_aitken(f, n0);
    n1 = min(2000, max(20, floor(Nmax / 2)));
    a2 = try_aitken(f, n1);

    if ~isnan(a1) && ~isnan(a2) && abs(a1 - a2) < 1e-6 * max(1, abs(a2))
        beta = a2;
    elseif ~isnan(a2)
        beta = a2;
    elseif ~isnan(a1)
        beta = a1;
    else
        % fallback: plain ratio at the largest available n
        if Nmax >= 2 && f(Nmax + 1) > 0
            beta = f(Nmax) / f(Nmax + 1);
        else
            beta = Inf;
        end
    end
end

function a = try_aitken(f, n0)
%TRY_AITKEN  Aitken extrapolation of the ratio sequence r(n) = f(n-1)/f(n)
%   at n0, 2*n0, 4*n0, assuming r(n) -> beta + O(1/n).
    Nmax = numel(f) - 1;
    n3 = 4 * n0;
    if n3 > Nmax || n3 < 1
        a = NaN;
        return;
    end
    r = @(n) ratio_at(f, n);
    r1 = r(n0); r2 = r(2*n0); r3 = r(n3);
    if any(~isfinite([r1 r2 r3]))
        a = NaN;
        return;
    end
    denom = r3 - 2*r2 + r1;
    if abs(denom) < 1e-13
        a = r3; % already converged / linear trend
        return;
    end
    extrap = r3 - (r3 - r2)^2 / denom;
    if isfinite(extrap) && extrap > 0
        a = extrap;
    else
        a = NaN;
    end
end

function r = ratio_at(f, n)
%RATIO_AT  f(n-1)/f(n) using 1-indexed storage (f(1) = f(0), f(n+1) = f(n)).
    if n < 1 || n + 1 > numel(f) || f(n + 1) <= 0
        r = NaN;
        return;
    end
    r = f(n) / f(n + 1);
end

function d = density_at_z(f, z, Nmax)
%DENSITY_AT_Z  rho(z) = z F'(z) / F(z), computed in log-space to avoid
%   overflow, matching the convention used elsewhere in this toolkit.
    logz = log(z);
    n = (0:Nmax)';
    fCol = f(:);
    validF = fCol > 0;

    logTermsF = n(validF) .* logz + log(fCol(validF));
    logF = logsumexp_local(logTermsF);

    nGE1 = validF & (n >= 1);
    logTermsFp = (n(nGE1) - 1) .* logz + log(fCol(nGE1)) + log(n(nGE1));
    logFp = logsumexp_local(logTermsFp);

    if ~isfinite(logF)
        d = NaN;
        return;
    end
    d = z * exp(logFp - logF);
end

function s = logsumexp_local(x)
%LOGSUMEXP_LOCAL  Numerically stable log(sum(exp(x))).
    if isempty(x)
        s = -Inf;
        return;
    end
    m = max(x);
    if ~isfinite(m)
        s = -Inf;
        return;
    end
    s = m + log(sum(exp(x - m)));
end
