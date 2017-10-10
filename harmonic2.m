function [value] = harmonic2(x)
    y = (1:x).^2;
    value = sum(1./y);
end
