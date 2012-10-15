function [positions, lengths, values] = boolean_domains(tests, is_circular)

  if (nargin == 1)
    is_circular = false;
  end

  if (isempty(tests))
    positions = NaN;
    lengths = NaN;
    values = NaN;
    
    return;
  end
  
  tests = tests(:);
  positions = find([1; diff(tests)]);
  lengths = diff([positions; length(tests)+1]);
  values = tests(positions);

  if (is_circular & (values(1) == values(end)) & numel(values) > 2)
    lengths(1) = lengths(1) + lengths(end);
    lengths(end) = lengths(1);
  end

  return;
end
