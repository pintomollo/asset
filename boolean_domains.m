function [positions, lengths, values] = boolean_domains(tests)

  tests = tests(:);
  positions = find([1; diff(tests)]);
  lengths = diff([positions; length(tests)]);
  values = tests(positions);

  return;
end
