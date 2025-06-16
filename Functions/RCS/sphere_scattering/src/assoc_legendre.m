function P = assoc_legendre(n, m, x)
% P = assoc_legendre(n, m, x) computes the associated Legendre function
% of degree n and order m. It is a wrapper function of MATLAB
% legendre(n, x) function.  Note that Matlab uses the definition of
% the associated Legendre function with Condon-Shortly phase factor.
% 
% n        an integer vector denoting the degree.
% m        an integer vector denoting the order. When m==0, we compute the
%          legendre's polynomial.  m must statisfy the condition m <= n.
% x        contain real values in [-1, 1] and must be a column vector.
%
% n and m cannot simultaneously be vectors.

    x = reshape(x, 1, length(x));
    
    if (length(n) > 1 && length(m) > 1)
        error(['Both n and m have length greater than ONE. Only one, either ' ...
               'the order or the degree can be a vector.']);
    end
    
    if (length(n) == 1 && length(m) >= 1)
        if n > max(m)
            a = legendre(n,x);
            P = a(m+1,:);
        else
            P = zeros(length(m), length(x));
            a = legendre(n,x);
            idx = find(m<=n);
            if (~isempty(idx))
                P = a(m(idx)+1,:);
            end
        end
    else
        if m > max(m)
            P = zeros(length(n), length(x));
        else
            for i = 1:length(n)
                P(i,:) = assoc_legendre(n(i), m, x);
            end
        end
    end
    
    
    