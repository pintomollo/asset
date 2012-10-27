function mystruct = func2mex(mystruct)

  switch class(mystruct)
    case 'handle'
      fname = func2str(mystruct);
      if (exist([fname '_mex'], 'file') == 3)
        mystruct = str2func([fname '_mex']);
      end
    case 'struct'
      fields = fieldnames(mystruct);
      nelems = numel(mystruct);

      for i=1:nelems
        for f=1:length(fields)
          mystruct(i).(fields{f}) = func2mex(mystruct(i).(fields{f}));
        end
      end
  end

  return;
end
