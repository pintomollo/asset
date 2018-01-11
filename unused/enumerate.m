function enum = enumerate(varargin)

  if (nargin == 1)
    enum = varargin{1}(:);
  else
    values = varargin{1};
    others = enumerate(varargin{2:end});

    nvalues = length(values);
    [nothers, ncols] = size(others);

    enum = NaN(nvalues*nothers, ncols+1);

    for i=1:nvalues
      indx = [((i-1)*nothers) + 1 : i*nothers];
      enum(indx, 1) = values(i);
      enum(indx, 2:end) = others;
    end
  end

  return;
end
