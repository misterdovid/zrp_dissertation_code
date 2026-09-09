function f = build_f_from_u(u_fn, Nmax)
%BUILD_F_FROM_U  Construct single-site weights f(0..Nmax) from a departure
%   rate function u(n).
%
%   f(0) = 1
%   f(n) = f(n-1) / u(n),   n = 1, ..., Nmax
%
%   u_fn must be a function handle accepting n = 1, 2, ..., Nmax (as a
%   scalar or vector) and returning positive values.

    f = zeros(1, Nmax + 1);
    f(1) = 1; % f(0) in 1-indexed MATLAB storage
    for n = 1:Nmax
        u = u_fn(n);
        if ~(u > 0) || ~isfinite(u)
            error('build_f_from_u:badRate', ...
                'u(%d) = %.6g is not a positive finite number.', n, u);
        end
        f(n + 1) = f(n) / u;
    end
end
