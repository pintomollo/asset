function spots = detect_spots(img, thresh, coef)

  if (nargin < 2)
    thresh = 0.05;
    coef = 3;
  elseif (nargin < 3)
    coef = 3;
  end
 
  rad_min = 1;
  atrous = imatrou(img);
  [h,w,nlayers] = size(atrous);

  %for i=1:nlayers
  %  figure;implot(atrous(:,:,i));
  %end

  nlayers = nlayers - 1;

  proj = prod(abs(atrous(:,:,2:end)),3);

  %figure;implot(proj);

  %proj = ones(h,w);
  %layers = zeros(h,w,nlayers);

  %for i=1:nlayers
  %  filtered = abs(atrous(:,:,i));
  %  mad_thresh = coef * mad(filtered(:),1) / 0.6745;
  %  filtered(filtered < mad_thresh) = 0;
  %  proj = proj .* filtered;
  %  layers(:,:,i) = filtered;
  %  %img = filtered .* img;
  %  %atrous(:,:,i) = (abs(img) > thresh); 
  %end

  %figure;implot(proj);

  %keyboard
 
  %oks = (proj > (thresh^nlayers));
  oks = (proj > coef * 1.253 * mad(proj(:)));
  spots = imfilter(double(oks), ones(3), 'symmetric');
  spots = (spots >= 7) & oks & imregionalmax(proj, 8);

  %figure;imshow(spots);

  [cand_x, cand_y] = find(spots);
  nspots = length(cand_x);

  %figure;imshow(img); hold on;
  opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');

  results = NaN(nspots, 5, nspots);

  for i=1:nspots
    [pos, sigma, ampl, bkg] = estimate_spot(img, [cand_x(i),cand_y(i)]);

    if (isempty(pos))
      break;
    end

    [best] = lsqnonlin(@fit_gaussian, [pos sigma ampl bkg], [], [], opts);

    if (i == 1)
      results(1,:,1) = best;
    else
      dists = hypot(results(:,1,1) - best(1), results(:,2,1) - best(2)) ;
      indx = find(dists < results(:,3,1) + best(3), 1);

      if (isempty(indx))
        results(find(isnan(results(:,1,1)),1),:,1) = best;
      else
        results(indx,:,find(isnan(results(indx,1,:)),1)) = best;
      end
    end
  end

  %keyboard

  results = mymean(results, 3);
  spots = results(~isnan(results(:,1)),:);
  spots = spots(spots(:,3) > rad_min,:);

  [junk, indx] = sort(spots(:,4),1,'descend');
  spots = spots(indx,:);

  %scatter(spots(:,2),spots(:,1),'r');

  %candidates = img(spots);
  %[candidates, indx] = sort(candidates,'descend');
  %cand_x = cand_x(indx);
  %cand_y = cand_y(indx);

  %waves = zeros(nspots, nlayers);
  %maxs = zeros(1,nspots);
  %for i=1:nspots
  %  waves(i,:) = layers(cand_y(i), cand_x(i), :);
  %  maxs(i) = find(waves(i,:)==max(waves(i,:),[],2),1);
  %end

  %gauss = GaussMask2D(spots(1,3), size(img), spots(1,1:2), 0, 1)*spots(1,4) + GaussMask2D(spots(2,3), size(img), spots(2,1:2), 0, 1)*spots(2,4) + mean(spots(1:2,5));
  %imshow(gauss)
  %atrous = mean(atrous(:,:,2:5),3);

  %keyboard

  return;

  function err = fit_gaussian(params)
    
    tmp_pos = params(1:2);
    tmp_sigma = params(3);
    tmp_ampl = params(4);
    tmp_bkg = params(5);

    pix_pos = round(tmp_pos);
    wsize = ceil(3*tmp_sigma);
    window = get_window(img, pix_pos, wsize);
    gauss = GaussMask2D(sigma, 2*wsize+1, pix_pos - tmp_pos);
    gauss = gauss*tmp_ampl + tmp_bkg;

    err = sum(sum(abs(gauss - window)));

    return;
  end
end

function [pos, sigma, ampl, bkg] = estimate_spot(img, estim_pos)

  %[h,w] = size(img);
  wsize = 50;
  %spots = NaN(2*wsize + 1);

  spots = get_window(img, estim_pos, wsize);

  %indx = [1:2*wsize+1];
  %indxx = indx - wsize - 1 + posx;
  %okx = (indxx > 0 & indxx <= h);
  %indxy = indx - wsize - 1 + posy;
  %oky = (indxy > 0 & indxy <= w);
  %spots(indx(okx), indx(oky)) = img(indxx(okx), indxy(oky));
  
  %mad_thresh = 2*mad(spots(:),1) / 0.6745;
  avg = mymean(spots(:));

  %projx = sum(spots,1);
  %projy = sum(spots,2).';

  %peakx = double(projx > mean(projx) + std(projx)) .* [-wsize:wsize];
  %peaky = double(projy > mean(projy) + std(projy)) .* [-wsize:wsize];
  
  bw = (spots > avg + 3*mad(spots(:)));
  bw = bwareaopen(bw, 5, 4);
  props = regionprops(bw, spots, 'Centroid', 'EquivDiameter', 'MaxIntensity', 'MinIntensity');

  %keyboard

  %if (~any(peakx) | ~any(peaky))
  if (length(props) == 0)
    pos = [];
    sigma = 0;
    ampl = 0;
    bkg = 0;

    return;
  elseif (length(props) == 1)
    indx = 1;
  else
    indx = 1;
    d = hypot(props(1).Centroid(1)-wsize-1,props(1).Centroid(2)-wsize-1);
    for i=2:length(props)
      tmp = hypot(props(i).Centroid(1)-wsize-1,props(i).Centroid(2)-wsize-1);
      if (tmp < d)
        d = tmp;
        indx = i;
      end
    end

  end

  pos = props(indx).Centroid;
  pos = pos(1, [2 1]);

  sigma = props(indx).EquivDiameter / 3;
  ampl = props(indx).MaxIntensity;
  bkg = props(indx).MinIntensity;
  ampl = ampl - bkg;
    
  gauss = GaussMask2D(sigma, 2*wsize+1, pos, 0, 1) * ampl + bkg;

  %figure;
  %subplot(121)
  %imshow(spots);
  %subplot(122)
  %imshow(gauss);
  %drawnow

  %tmp1 = find(peakx > 0, 1, 'first');
  %tmp2 = find(peakx < 0, 1, 'last');
  %pos(1) = tmp1 + tmp2 / 2;
  %sigma = tmp1 - tmp2 / 2;

  %tmp1 = find(peaky > 0, 1, 'first');
  %tmp2 = find(peaky < 0, 1, 'last');
  %pos(2) = tmp1 + tmp2 / 2;
  %sigma = sigma +  ((tmp1 - tmp2) / 2) / 2;

  %[x, y] = find(spots > avg + mad_thresh);
  %[pos, sigma] = mymean([x y], [], 1);
  %sigma = mean(sigma.');

  pos = pos - wsize - 1 + estim_pos;

  return;
end

function window = get_window(img, pos, wsize)

  [h,w] = size(img);
  window = NaN(2*wsize + 1);

  indx = [1:2*wsize+1];
  indxx = indx - wsize - 1 + pos(1);
  okx = (indxx > 0 & indxx <= h);
  indxy = indx - wsize - 1 + pos(2);
  oky = (indxy > 0 & indxy <= w);

  window(indx(okx), indx(oky)) = img(indxx(okx), indxy(oky));
  
  return;
end
