function y = norm_nan(x)
    id = isnan(x);
    y = norm(x(~id));
end