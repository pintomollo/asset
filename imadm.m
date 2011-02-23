function [edge, direct] = imadm(img, thresh, single)
  % IMADM implements the Absolute Difference Mask edge detector
  %     [EDGE, DIRECT] = IMADM(IMG) computes the strength and the
  %     direction of the edges (EDGES and DIRECT) of the image. The
  %     strength ranges from 0 to 1, the direction follows the convention :
  %       1 = Negative diagonal (antidiagonal)
  %       2 = Vertical
  %       3 = Positive diagonal
  %       4 = Horizontal
  %
  %     [...] = IMDAM(IMG, THRESH) imposes a threshold on the
  %     value of the edges (standard is 0).
  %
  %     EDGE = IMDAM(...) returns only the edge intensity map.
  %
  % References:
  % [1] F.M. Alzahrani, T. Chen, "A real-time edge detector: algorithm and VLSI architecture", Real-Time Imaging 3 (1997) 363-378. 
  % [2] S.-C. Zhang and Z.-Q. Liu, "A robust, real-time ellipse detector," Pattern Recognition, vol. 38, no. 2, (2005) 273-287.

  if(nargin < 2)
    thresh = 0;
    single = true;
  elseif (nargin < 3)
    single = true;
  end

  if (isempty(thresh))
    thresh = 0;
  end

  % Circular Gaussian filtre as defined in [2]
  avgfilt = [0.25 0.5  0.5 0.5  0.25; ...
             0.5  0.75 1   0.75 0.5; ...
             0.5  1    2   1    0.5; ...
             0.5  0.75 1   0.75 0.5; ...
             0.25 0.5  0.5 0.5  0.25] / 16;

  img = imfilter(img, avgfilt, 'symmetric');

  % The filters of the ADM as defined in [1]
  hfilt = [-1 -1 0 1 1];
  vfilt = -hfilt';
  pdfilt = -diag(hfilt);
  ndfilt = fliplr(pdfilt);
  
  % Compute the edges along the 4 different directions
  himg = imfilter(img, hfilt, 'symmetric');
  vimg = imfilter(img, vfilt, 'symmetric');
  pdimg = imfilter(img, pdfilt, 'symmetric');
  ndimg = imfilter(img, ndfilt, 'symmetric');

  % Create a structure to find easily the maxima and minima of every pixel
  img = abs(cat(3, ndimg, vimg, pdimg, himg));

  edge = max(img, [], 3) / 2;
  [tmp, direct] = min(img, [], 3);

  % Create a framework necessary to compare each cell to its 8 neighbors
  frame = ones(size(edge) + 2);
  frame(2:end-1, 2:end-1) = edge;

  % Reduce each edge to a single pixel-wide edge
  singleedge = zeros(size(edge));
  for i=1:4
    indx = (direct==i);

    % The flowchart as defined in [1], keeps only the local maximas
    switch i
      case 1,
        compar = (edge >= frame(1:end-2, 1:end-2) & edge > frame(3:end, 3:end));
      case 2,
        compar = (edge >= frame(2:end-1, 1:end-2) & edge > frame(2:end-1, 3:end));
      case 3,
        compar = (edge >= frame(1:end-2, 3:end) & edge > frame(3:end, 1:end-2));
      case 4,
        compar = (edge >= frame(1:end-2, 2:end-1) & edge > frame(3:end, 2:end-1));
    end 

    indx = (indx & compar);
    singleedge(indx) = edge(indx);
  end

  if (single)
    edge = singleedge;
  end

  minimg = min(min(edge));
  maximg = max(max(edge));

  edge = (edge - minimg) / (maximg - minimg);

  % Apply the threshold
  edge(edge < thresh) = 0;

  if (nargout == 2)
    % Compute the edge angle based on the intensities
    direct = -atan2(himg, vimg);
    tmpdirect = -(atan2(ndimg, pdimg) - pi/4);
    probs = (direct > pi/2 & tmpdirect < -pi/2);
    tmpdirect(probs) = tmpdirect(probs) + 2*pi;
    probs = (tmpdirect > pi/2 & direct < -pi/2);
    direct(probs) = direct(probs) + 2*pi;
    direct = (direct + tmpdirect) ./ 2;
    direct(direct < 0) = direct(direct < 0) + 2*pi;

    direct(edge == 0) = 0;
  end

  return;
end
