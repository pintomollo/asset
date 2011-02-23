function struct_size(mystruct, indent, mult, total)

  if (nargin < 2)
    myfield = mystruct;
    tmp = whos('myfield');
    total = tmp.bytes;

    indent = '|_';
    mult = 1;

    display(['TOTAL: ' num2str(total)]);
  end

  if (isstruct(mystruct))

    names = fieldnames(mystruct);

    for i=1:length(names)
      myfield = mystruct.(names{i});
      tmp = whos('myfield');
      nelems = prod(tmp.size);
      display([indent names{i} ':' num2str(100 * (tmp.bytes * mult) / total) '%'])
      if (isstruct(myfield))
        struct_size(mystruct.(names{i})(1), [' ' indent], mult * nelems, total);
      end

    end
  else
    myfield = mystruct;
    tmp = whos('myfield');
    sizes = tmp.bytes;
  end

  return;
end
