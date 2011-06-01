function [datas, bkg, cyto, theta] = gather_quantification(mymovie, opts)

  npts = 500;

  nframes = size_data(mymovie.data);  
  theta = [-0.5:1/npts:0.5].';
  theta = theta(1:end-1);

  datas = NaN(nframes, npts);
  bkg = NaN(nframes, 7);

  type = 'linear';

  switch type
    case 'elliptic'
      theta = theta * 2*pi;
      range = [-pi pi];
    case 'linear'
      range = [-0.5 0.5];
  end

  for i=1:nframes
    warper = mymovie.data.warpers(i);

    if (isfield(mymovie.data.quantification, 'carth') & ~isempty(mymovie.data.quantification(i).carth))
      cortex = mymovie.data.quantification(i).carth;
    else
      cortex = mymovie.data.cortex(i).carth;
    end
    cortex = carth2normalized(cortex, warper, opts);
    switch type
      case 'elliptic'
        pts = carth2elliptic(cortex, warper.reference.centers, warper.reference.axes_length, warper.reference.orientations);
        pts = pts(:, 1);
      case 'linear'
        %pts = carth2linear(cortex, mymovie.markers.ruffles(i).carth, mymovie.markers.ruffles(i).paths)
        pts = carth2linear(cortex);
    end
    %full_theta = mymovie.data.warpers(i).warp(:,1);
    intens = mymovie.data.quantification(i).cortex;

    new_intens = interp_elliptic([pts, intens], theta, range);
    datas(i, :) = new_intens(:, 2);
    %bkg(i,:) = mymovie.data.quantification(i).bkg;
  end

  if (isfield(mymovie.data, 'cytokinesis'))
    cyto = mymovie.data.cytokinesis;
  else
    cyto = 0;
  end

  return;
end
