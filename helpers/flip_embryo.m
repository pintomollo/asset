function mymovie = flip_embryo(mymovie, opts)

  type = opts.segmentation_type;
  orient = mymean(mymovie.(type).orientations) + pi;
  orient = orient - 2*pi*(orient > 2*pi);

  fields = fieldnames(mymovie);
  for f = 1:length(fields)
    if (~isempty(mymovie.(fields{f})) && isfield(mymovie.(fields{f}), 'orientations') && ~isempty(mymovie.(fields{f}).orientations))

      mymovie.(fields{f}) = align_orientations(mymovie.(fields{f}), orient);
      mymovie.(fields{f}).inverted = true;
    end
  end

  if (isfield(mymovie, 'metadata') & ~isempty(mymovie.metadata) & isfield(mymovie.metadata, 'orientation_3d') & ~isempty(mymovie.metadata.orientation_3d))
    mymovie.metadata.orientation_3d = align_orientations(mymovie.metadata.orientation_3d, orient);
  end

  return;
end
