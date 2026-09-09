function p = pGCE_from_z(f, z, Nmax)
%PGCE_FROM_Z  Grand-canonical single-site marginal p(n; z) = z^n f(n) / F(z)
%   for n = 0..Nmax, computed in log-space for numerical stability.

    n = (0:Nmax)';
    fCol = f(:);
    if numel(fCol) < Nmax + 1
        error('pGCE_from_z:tooShort', 'f must have at least Nmax+1 entries.');
    end
    fCol = fCol(1:Nmax + 1);

    logz = log(z);
    validF = fCol > 0;
    logTerms = n(validF) .* logz + log(fCol(validF));
    m = max(logTerms);
    logF = m + log(sum(exp(logTerms - m)));

    p = zeros(Nmax + 1, 1);
    p(validF) = exp(n(validF) .* logz + log(fCol(validF)) - logF);
end
