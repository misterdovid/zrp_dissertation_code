function [coeffs, logScale] = poly_pow_trunc_scaled(a, L, maxDeg)
%POLY_POW_TRUNC_SCALED  Compute F(z)^L truncated to degree maxDeg, where
%   F(z) has coefficients a(1), a(2), ... (a(1) = coefficient of z^0),
%   via binary exponentiation (repeated squaring).
%
%   Returns coeffs (1 x (maxDeg+1), coefficients of z^0 .. z^maxDeg) and a
%   scalar logScale such that the TRUE coefficients are
%       true_coeffs = coeffs * exp(logScale)
%
%   Rescaling: Z_{L,N} = [z^N] F(z)^L grows roughly like c^L for fixed
%   density N = rho*L (it's an extensive thermodynamic quantity), which
%   overflows double precision (~1e308) once L reaches a few thousand at
%   typical densities. We periodically renormalize the coefficient array
%   by its own maximum whenever that maximum's magnitude approaches the
%   overflow threshold, and track the accumulated log-correction
%   separately. This keeps all arithmetic in plain, fast, unscaled
%   double precision while still reaching very large L.
%
%   Deliberately uses plain O(N^2) convolution (MATLAB's CONV), not FFT:
%   the intermediate coefficients here can span 10+ orders of magnitude
%   within a single array, and double-precision FFT convolution loses
%   catastrophic relative precision on inputs like that (verified during
%   development). Plain convolution has no such failure mode.

    aRow = a(:)'; % ensure row vector
    if numel(aRow) < maxDeg + 1
        aRow = [aRow, zeros(1, maxDeg + 1 - numel(aRow))];
    else
        aRow = aRow(1:maxDeg + 1);
    end

    result = zeros(1, maxDeg + 1);
    result(1) = 1; % z^0 coefficient of F^0 = 1
    resultLogScale = 0;

    base = aRow;
    baseLogScale = 0;

    e = L;
    while e > 0
        if bitand(e, 1) == 1
            result = poly_mul_trunc(result, base, maxDeg);
            resultLogScale = resultLogScale + baseLogScale;
            [result, adj] = renormalize(result);
            resultLogScale = resultLogScale + adj;
        end
        e = bitshift(e, -1);
        if e > 0
            base = poly_mul_trunc(base, base, maxDeg);
            baseLogScale = 2 * baseLogScale;
            [base, adj] = renormalize(base);
            baseLogScale = baseLogScale + adj;
        end
    end

    coeffs = result;
    logScale = resultLogScale;
end

function c = poly_mul_trunc(a, b, maxDeg)
%POLY_MUL_TRUNC  Truncated polynomial multiplication via plain convolution.
    full = conv(a, b);
    len = min(numel(full), maxDeg + 1);
    c = zeros(1, maxDeg + 1);
    c(1:len) = full(1:len);
end

function [arrOut, logm] = renormalize(arr)
%RENORMALIZE  If the array's maximum magnitude is large enough that
%   continued growth risks double-precision overflow, rescale the array
%   so its max is O(1) and return the log of the factor divided out. If
%   there's still plenty of headroom, return the array unchanged and
%   logm = 0. Multiply arrOut by exp(logm) to recover the true values.
    m = max(arr);
    if m == 0 || ~isfinite(m)
        arrOut = arr;
        logm = 0;
        return;
    end
    logAbsM = log(m);
    if abs(logAbsM) < 200
        arrOut = arr;
        logm = 0;
        return;
    end
    arrOut = arr * exp(-logAbsM);
    logm = logAbsM;
end
