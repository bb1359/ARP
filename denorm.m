function params = denorm(x,exponent,offset)
  params = real(10.^(x.*exponent + offset));
  %params = (log10(x) - offset) ./ exponent;
end
