function channels = get_channels(mymovie)

  fields = fieldnames(mymovie);
  nfields = length(fields);
  channels = cell(nfields, 1);

  for i=1:nfields
    field = fields{i};

    if (~isempty(mymovie.(field)) & isfield(mymovie.(field), 'fname') & exist(mymovie.(field).fname, 'file'))
      channels{i} = field;
    end
  end

  channels = channels(~cellfun('isempty', channels));

  return;
end
