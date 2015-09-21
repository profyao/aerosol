function z = corr_nan(x,y)
        id = ~isnan(x) & ~isnan(y);
        z = corr(x(id),y(id));
end