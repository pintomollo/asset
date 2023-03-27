function display_graphs(fname, field)

  loading = false;
  if (nargin == 0)
    field = 'dic';
    loading = true;
  elseif (nargin == 1)
    if (exist(fname, 'file') | exist([fname '.mat'], 'file'))
      field = 'dic';
    else
      field = fname;
      loading = true;
    end
  end

  if (loading && exist('results_segment.mat', 'file'))
    load('results_segment');
    opts = get_struct('RECOS',1);
    opts.recompute = true;
  else
    load(fname);
    opts = get_struct('RECOS',1);
    switch field
      case 'dic'
        [imgsize nframes] = size_data(mymovie.dic);

        opts.segmentation_parameters = set_image_size(opts.segmentation_parameters, imgsize);

        opts.measure_performances = true;
        opts.analysed_fields = {'estim', 'carth'};
        opts.recompute = true;

        mymovie = segment_movie(mymovie, opts);
        mymovie2 = tracing_error(trackings, mymovie, field, opts);

        opts = load_parameters('optimized_params', opts);
        opts.measure_performances = false;
        opts.analysed_fields = {'carth'};
        opts.segmentation_parameters = set_image_size(opts.segmentation_parameters, imgsize);

        mymovie = segment_movie(mymovie, opts);
        mymovie = tracing_error(trackings, mymovie, field, opts);
      
        save('results_segment','mymovie*','trackings');

      case 'markers'
        [imgsize nframes] = size_data(mymovie.cortex);
        opts.segmentation_type = 'dic';
        opts = load_parameters('optimized_params', opts);
        opts.recompute = true;

        mymovie = segment_movie(mymovie, opts);
        mymovie.markers = mymovie.dic;
        mymovie2 = tracing_error(trackings, mymovie, field, opts);

        opts.segmentation_type = 'markers';
        mymovie = segment_movie(mymovie, opts);
        %mymovie = correct_dic_shift(mymovie, opts.segmentation_parameters.correction, opts);
        mymovie3 = tracing_error(trackings, mymovie, field, opts);

        mymovie.markers.cortex = [];
        opts.recompute = false;
        opts.segmentation_type = 'markers';
        mymovie = segment_movie(mymovie, opts);
        mymovie4 = tracing_error(trackings, mymovie, field, opts);

        mymovie.markers.cortex = [];
        opts = get_struct('RECOS',1);
        opts.segmentation_type = 'markers';
        opts.measure_performances = true;
        opts.analysed_fields = {'estim', 'carth'};
        mymovie = segment_movie(mymovie, opts);
        mymovie = tracing_error(trackings, mymovie, field, opts);

        save('results_segment','mymovie*','trackings');

      otherwise
        disp(['Unknown type: ' field]);
        return;
    end
  end

  %compare_manuals(mymovie, trackings, field, opts);

  switch field
    case 'dic'
      orig_errors = mymovie2.(field).errors;
      estim_errors = orig_errors(:,:,1,:);
      orig_errors = orig_errors(:,:,2,:);
      optim_errors = mymovie.(field).errors;
      track_errors = trackings.(field).errors;

      [egg, cortex, index_egg, index_cortex] = compute_errors(estim_errors);
      [egg, cortex, index_egg, index_cortex] = compute_errors(orig_errors, egg, cortex, index_egg, index_cortex);
      [egg, cortex, index_egg, index_cortex] = compute_errors(optim_errors, egg, cortex, index_egg, index_cortex);
      [egg, cortex, index_egg, index_cortex] = compute_errors(track_errors, egg, cortex, index_egg, index_cortex);

      names = {'Initial', 'ASSET', 'ASSET + ML', 'Manual'};

      [H,p] = myttest(egg, index_egg)
      [H,p] = myttest(cortex, index_cortex)
      %[H,p] = myttest([egg;cortex], [index_egg; index_cortex + max(index_egg)])

      [errs, stds] = mymean(egg, [], index_egg);
      [errs(:,2), stds(:,2)] = mymean(cortex, [], index_cortex);
      figure;
      barweb(errs.' * 100, stds.' * 100, [], {'Eggshell', 'Cortex'}, '', 'Segmentation type', 'Average distance to the reference (% of embryo''s radius)',[],'y',[],1);

    case 'markers'
      orig_errors = mymovie.(field).errors;

      estim_errors = orig_errors(:,:,1,:);
      orig_errors = orig_errors(:,:,2,:);
      optim_errors = mymovie4.(field).errors;
      corr_errors = mymovie3.(field).errors;
      dic_errors = mymovie2.(field).errors;
      track_errors = trackings.(field).errors;

      [egg, cortex, index_egg, index_cortex] = compute_errors(estim_errors);
      [egg, cortex, index_egg, index_cortex] = compute_errors(dic_errors, egg, cortex, index_egg, index_cortex);
      [egg, cortex, index_egg, index_cortex] = compute_errors(corr_errors, egg, cortex, index_egg, index_cortex);
      [egg, cortex, index_egg, index_cortex] = compute_errors(orig_errors, egg, cortex, index_egg, index_cortex);
      [egg, cortex, index_egg, index_cortex] = compute_errors(optim_errors, egg, cortex, index_egg, index_cortex);
      [egg, cortex, index_egg, index_cortex] = compute_errors(track_errors, egg, cortex, index_egg, index_cortex);

      names = {'Initial', 'DIC', 'Correction', 'mCherry', 'mCherry + ML', 'Manual'};

      [H,p] = myttest(cortex, index_cortex)
      %[H,p] = myttest([egg;cortex], [index_egg; index_cortex + max(index_egg)])

      [errs, stds] = mymean(cortex, [], index_cortex);
      figure;
      barweb(errs.' * 100, stds.' * 100, [], '', '', 'Segmentation type', 'Error (% of ER)',[],'y',[],1);

  end

  display_angular;

  return;
end

function compare_manuals(mymovie, trackings, field, opts)

  optim_errors = mymovie.(field).errors;
  track_errors = trackings.(field).errors;

  [egg, cortex, index_egg, index_cortex] = compute_errors(optim_errors);
  [egg, cortex, index_egg, index_cortex] = compute_errors(track_errors, egg, cortex, index_egg, index_cortex);

  nmanuals = length(trackings.(field).child);
  for i=1:nmanuals
    manual_error(i).err = spline_error(trackings.(field).child(i), mymovie.(field), opts.nbins, opts);
    splines = [];
    for j=1:nmanuals
      if (j~=i)
        splines = cat(3, splines, trackings.(field).child(j).mean);
      end
    end
    inter_error(i).err = spline_error(trackings.(field).child(i), splines, opts.nbins, opts);
    [egg, cortex, index_egg, index_cortex] = compute_errors(manual_error(i).err, egg, cortex, index_egg, index_cortex);
    [egg, cortex, index_egg, index_cortex] = compute_errors(inter_error(i).err, egg, cortex, index_egg, index_cortex);
  end

  [errs, stds] = mymean(egg, [], index_egg);
  errs = reshape(errs, 2, []);
  stds = reshape(stds, 2, []);
  [tmp_errs, tmp_stds] = mymean(cortex, [], index_cortex);
  tmp_errs = reshape(tmp_errs, 2, []);
  tmp_stds = reshape(tmp_stds, 2, []);
  errs = cat(1, errs, tmp_errs);
  stds = cat(1, stds, tmp_stds);
  figure;
  barweb(errs.' * 100, stds.' * 100, [], {'Reference', 'Manual #A', 'Manual #B', 'Manual #C', 'Manual #D' }, '', 'Segmentation type', 'Average distance to the reference (% of embryo''s radius)',[],'y',[],1);

  [H,p] = myttest([egg;cortex], [index_egg; index_cortex + max(index_egg)])

  return;
end

function [eggs, cors, indxe, indxc] = compute_errors(errors, eggs, cors, indxe, indxc)

  egg = errors(1,:,:,:);
  egg = egg(:);

  cor = errors(2,:,:,:);
  cor = cor(:);

  if (nargin > 1)
    eggs = [eggs;egg];
    cors = [cors; cor];

    indxe = [indxe; ones(size(egg))*(max(indxe) + 1)];
    indxc = [indxc; ones(size(cor))*(max(indxc) + 1)];
  else
    eggs = egg;
    cors = cor;

    indxe = ones(size(egg));
    indxc = ones(size(cor));
  end


  return;
end
