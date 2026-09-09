function [p, meanUCE] = pCE_repeated_squaring(f, L, N)
%PCE_REPEATED_SQUARING  Canonical single-site marginal p_CE(0..N) for a
%   homogeneous ZRP with L sites and N particles, computed via repeated
%   squaring of the truncated generating polynomial (fast, no FFT; see
%   poly_pow_trunc_scaled.m for why FFT is avoided here).
%
%   f        : row or column vector of single-site weights f(0), f(1), ...
%              (length must be at least N+1)
%   L, N     : positive integers, number of sites and particles
%
%   Returns:
%   p        : column vector, p(n+1) = p_CE(n) for n = 0..N
%   meanUCE  : <u(n)>_CE = Z_{L,N-1} / Z_{L,N}  (Evans-Hanney Eq. 11);
%              this should converge to the exact GCE fugacity z(rho) as
%              L -> infinity (Eq. 37) -- a useful numerical check of
%              finite-L convergence.

    [ZL, logScaleL] = poly_pow_trunc_scaled(f, L, N);
    if ~(ZL(N + 1) > 0) || ~isfinite(ZL(N + 1))
        error('pCE_repeated_squaring:overflow', ...
            ['Z_{L,N} coefficient is not a usable positive finite number ' ...
             'even after rescaling. This should not happen for typical ' ...
             'inputs; check f(n) for validity.']);
    end
    [ZLm1, logScaleLm1] = poly_pow_trunc_scaled(f, L - 1, N);

    scaleDiff = logScaleLm1 - logScaleL;
    logZLN = log(ZL(N + 1));

    p = zeros(N + 1, 1);
    fRow = f(:)';
    for n = 0:N
        idx = N - n; % 0-indexed target power of the L-1 polynomial
        idxMat = idx + 1; % 1-indexed into ZLm1
        if idx < 0 || idxMat > numel(ZLm1) || ~(ZLm1(idxMat) >= 0) || ~(fRow(n + 1) > 0)
            p(n + 1) = 0;
            continue;
        end
        logRatio = log(ZLm1(idxMat)) - logZLN + scaleDiff;
        p(n + 1) = exp(log(fRow(n + 1)) + logRatio);
    end

    meanUCE = NaN;
    if N >= 1 && ZL(N) >= 0 % ZL(N) is the 1-indexed slot for exponent N-1
        meanUCE = ZL(N) / ZL(N + 1);
    end
end
