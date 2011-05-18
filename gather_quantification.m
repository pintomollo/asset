function [datas, theta, cyto] = gather_quantification(mymovie, opts)

  npts = 500;

  nframes = size_data(mymovie.data);  
  theta = [-pi:2*pi/npts:pi].';
  theta = theta(1:end-1);

  datas = NaN(nframes, npts);

  for i=1:nframes
    keyboard

    warper = mymovie.data.warpers(i);

    cortex = mymovie.data.cortex(i).carth;
    cortex = carth2normalized(cortex, warper, opts);
    ell_pts = carth2elliptic(cortex, warper.reference.centers, warper.reference.axes_length, warper.reference.orientations);
    %full_theta = mymovie.data.warpers(i).warp(:,1);
    intens = mymovie.data.quantification(i).cortex;

    new_intens = interp_elliptic([ell_pts(:,1), intens], theta);
    datas(i, :) = new_intens(:, 2);
  end

  if (isfield(mymovie.data, 'cytokinesis'))
    cyto = mymovie.data.cytokinesis;
  else
    cyto = 0;
  end

  return;
end
