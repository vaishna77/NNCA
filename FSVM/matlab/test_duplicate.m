function [label] = test_duplicate(x)
    if (dpendent_fun(x))
        label = 1;
    else
        label = -1;
    end
end

