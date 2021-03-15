% [ni, nj] = squaredFactors(number, landscape)
% Return two closest divisors of a number.

% 2021-03-15. Leonardo Molina.
% 2021-03-15. Last modified.
function [ni, nj] = squaredFactors(number, landscape)
    landscape = nargin == 1 || landscape;
    number = round(number);
    x = 1:number;
    rr = x(~(rem(number, x)));
    cc = number ./ rr;
    [~, k] = mink(abs(rr - cc), 1);
    if landscape == (cc(k) < rr(k))
        ni = cc(k);
        nj = rr(k);
    else
        ni = rr(k);
        nj = cc(k);
    end
end