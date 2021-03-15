function out = clip(in, min, max)

    if in < min
        out = min;
    elseif in > max
        out = max;
    else
        out = in;
    end

end