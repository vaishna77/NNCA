function [out] = kernel(xs, x,Qchoice)
    r = (xs-x).^2;
    if(Qchoice == 0)
        out = exp(-sqrt(sum(r)));
    else
        out = exp(-sum(r));
    end
end
