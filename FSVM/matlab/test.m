function [out,out2] = test(xs, xs_in, ys, ys_in, alpha, alpha_in, b, x, Qchoice)
    [xs_N,temp] = size(xs);
    [xs_in_N,temp] = size(xs_in);
    out = 0;
    for i=1:xs_N
        out = out + alpha(i)*ys(i)*kernel(xs(i,:),x, Qchoice);
    end
    for i=1:xs_in_N
        out = out + alpha_in(i)*ys_in(i)*kernel(xs_in(i,:),x, Qchoice);
    end
    out = out + b;
    if (out>=0)
        out2 = 1;
    else
        out2 = -1;
    end
%     return out;
end

