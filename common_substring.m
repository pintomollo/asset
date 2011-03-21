function [new_str, end_str] = common_substring(str1, str2, replacement_str)

  new_str = '';
  if (isempty(str1) | isempty(str2))
    return;
  end

  nstr = min(length(str1), length(str2));

  sindx = NaN;
  eindx = NaN;
  for i=1:nstr
    if (isnan(sindx) & str1(i) ~= str2(i))
      sindx = i-1;
    end
    if (isnan(eindx) & str1(end-i+1) ~= str2(end-i+1))
      eindx = i-2;
    end
  end

  if (isnan(sindx))
    sindx = nstr;
  end
  if (isnan(eindx))
    eindx = nstr-1;
  end
  if (eindx < 0)
    eindx = 0;
  end
    
  sep1 = '';
  sep2 = '';
  if (sindx > 0 & str1(sindx) ~= '_')
    sep1 = '_';
  end
  if (eindx < nstr & str1(end-eindx) ~= '_')
    sep2 = '_';
  end
      
  if (nargin < 3)
    new_str = str1(1:sindx);
  else
    new_str = [str1(1:sindx) sep1 replacement_str sep2 str1(end-eindx:end)];
  end

  if (nargout > 1)
    end_str = str1(end-eindx:end);
  end

  return;
end
