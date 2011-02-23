function y = logistic(x, center, width, slope, background, amplitude, polar)

  if (nargin < 7)
    polar = true;
  end
  if (nargin == 2)
    width = center(2);
    slope = center(3);
    background = center(4);
    amplitude = center(5);
    center = center(1);
  end

  if (polar)
    if (center > pi)
      below = (x < center & x > center - pi);
      around = (x < center - pi);
    else
      below = ~(x > center & x < center + pi);
      around = (x > center + pi);
    end
  end

  y = NaN(size(x));
  y(below) = amplitude ./ (1 + exp(-slope*(x(below) - (center - width)))) + background;
  y(below & around) = amplitude ./ (1 + exp(-slope*(x(below & around) - (center + 2*pi - width)))) + background;
  y(~below) = amplitude ./ (1 + exp(+slope*(x(~below) - (center + width)))) + background;
  y(~below & around) = amplitude ./ (1 + exp(+slope*(x(~below & around) - (center + 2*pi + width)))) + background;

  return;
end
