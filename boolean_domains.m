function [positions, lengths, values] = boolean_domains(tests)

  tests = tests(:);
  positions = find([1; diff(tests)]);
  lengths = diff([positions; length(tests)+1]);
  values = tests(positions);

  return;
end
