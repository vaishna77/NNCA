function [out] = label_func(x)
    r = sum(x.^2);
    value = exp(-sqrt(r));
    if (value < x(1)+0.6)
        out = 1;
    else
        out = -1;
    end
end

