function [starts,ends,angles] = detect_ellipse(img, precision, st, ed, ang)

  %close all

  [edge, direct] = imadm_mex(img);

  precision = precision * 2 * pi;

   %nconv = 5;
  %hconvex = ones(nconv);
  %hconvex(:,ceil(nconv/2)+1:end) = -1;
  %vconvex = hconvex';
  %pdconvex = ones(nconv);
  %pdconvex = pdconvex-2*tril(ones(nconv),-1);
  %ndconvex = fliplr(pdconvex);

  %himg = sign(imfilter(edge, hconvex, 'symmetric'));
  %vimg = sign(imfilter(edge, vconvex, 'symmetric'));
  %pdimg = sign(imfilter(edge, pdconvex, 'symmetric'));
  %ndimg = sign(imfilter(edge, ndconvex, 'symmetric'));

  %figure;h=pcolor(himg)
  %set(h,'EdgeColor','none')

  thresh = graythresh(edge);
  binedge = double(edge>thresh);
 
  neighfilt = [0 1 0; 1 0 1; 0 1 0];
  connectivity = imfilter(binedge, neighfilt, 'symmetric');

  binedge(connectivity==0) = 0;
  binedge(connectivity>=3) = 0;
  edge(binedge==0) = 0;
  direct(binedge==0) = 0;

  %himg = sign(imfilter(edge, hconvex, 'symmetric'));
  %vimg = sign(imfilter(edge, vconvex, 'symmetric'));
  %pdimg = sign(imfilter(edge, pdconvex, 'symmetric'));
  %ndimg = sign(imfilter(edge, ndconvex, 'symmetric'));

  figure;implot(direct)
  convex = compute_convexity(binedge, direct, 5);
  figure;implot(convex)
  convex = direct + sign(convex)*pi/2;
  figure;implot(convex)
  return
  %figure;implot(convex)

  %figure;h=pcolor(convex);
  %set(h,'EdgeColor','none')
  %figure;
  %implot(floatdirect)
  %figure;h=pcolor(floatdirect);
  %set(h,'EdgeColor','none')
  %figure;imshow(edge)
  %convex(direct==1) = 

  %figure;h=pcolor(abs(himg-himg2))
  %set(h,'EdgeColor','none')

  [starts, ends, angles] = get_point_pairs(logical(binedge), direct, convex, precision);

  scatter(starts(1,:),starts(2,:),'r')
  return

  [pts] = get_421_points(starts, ends, angles, precision, binedge);

  implot(binedge);
  hold on;
  scatter(pts(1,:),pts(2,:),'r')
  scatter(pts(4,:),pts(5,:),'g')
  scatter(pts(7,:),pts(8,:),'b')
  scatter(pts(10,:),pts(11,:),'k')

  return;
end

function [pts] = get_421_points(starts, ends, angles, precision, binedge)

  neigh_size = 5; 
  pts = zeros(12,0);

  for i=1:length(angles)

    new_indx = find_matching_pair(starts, angles, i, neigh_size, precision);

    if(~isempty(new_indx))
      new_pts = get_421_data(starts, ends, i, new_indx);
      pts = [pts new_pts];
    end
  end

  return;
end

function [candidates] = find_matching_pair(pts, angles, indx, maxdist, precision)

  candidates = zeros(1,0);

  pt = pts(:,indx);
  angl = angles(1,indx)
  indexes = [1:length(angles)];
  
  dist = sqrt(((pts(1,:) - pt(1)).^2) + ((pts(2,:) - pt(2)).^2));
  isvalid = ((dist<maxdist) & (dist~=0));

  pts = pts(:,isvalid)
  angles = angles(1,isvalid)
  indexes = indexes(isvalid);

  scatter(pts(1,:),pts(2,:),'y');

  line = get_line(pt(1:2,1), angl, -maxdist, maxdist);
  common = intersect(line',pts(1:2,:)','rows');

  isvalid = check_tangents(angl, angles, precision)

  if(isempty(common) & sum(isvalid)>0)
    pts = pts(:,isvalid);

    angles = angles(1,isvalid);
    indexes = indexes(isvalid);

    isvalid = logical(zeros(1,length(angles)));
    iter = [1:length(isvalid)];
    for i=iter
      isconvex = check_convexity(pt(1:2,1), pts(1:2,i), pt(3,1), pts(3,i), precision);
      line = get_line(pts(1:2,i), angles(i), -maxdist, maxdist);
      common = intersect(line',pts(1:2,(iter~=i))','rows');
      isvalid(i) = (isconvex & isempty(common));
    end

    if(any(isvalid))
      scatter(pts(1,isvalid),pts(2,isvalid),'g');
      candidates = indexes(isvalid);
    end
  end

  return;
end

function new_pts = get_421_data(starts, ends, indx, new_indx)

  new_pts = zeros(12,length(new_indx));

  line = [starts(1:3,indx);ends(1:3,indx)];
  for i=1:length(new_indx)
    new_pts(:,i) = [line; starts(1:3,new_indx(i)); ends(1:3,new_indx(i))];
  end

  return
end

function [starts, ends, angles] = get_point_pairs(edge, tangent, convex, precision)

  starts = zeros(4,0);
  ends = zeros(4,0);
  angles = zeros(1,0);

  xpos = find(any(edge,1));
  ypos = find(any(edge,2)');

  xaxe = (max(xpos) - min(xpos));
  yaxe = (max(ypos) - min(ypos));

  mindist = min(xaxe, yaxe)/4;
  maxdist = max(xaxe, yaxe);

  implot(flipud(tangent))
  hold on
  [new_starts, new_ends, new_angles] = get_pairs([584;245], convex(245,584)-pi/4, mindist, maxdist, edge, tangent, convex, precision);
  kk

  for x=xpos
    for y=ypos
      if(edge(y,x))

        [new_starts, new_ends, new_angles] = get_pairs([x;y], convex(y,x)+pi/4, mindist, maxdist, edge, tangent, convex, precision);
        if(~isempty(new_starts))
          starts = [starts new_starts];
          ends = [ends new_ends];
          angles = [angles new_angles];
        end

        [new_starts, new_ends, new_angles] = get_pairs([x;y], convex(y,x)-pi/4, mindist, maxdist, edge, tangent, convex, precision);

        if(~isempty(new_starts))
          starts = [starts new_starts];
          ends = [ends new_ends];
          angles = [angles new_angles];
        end

      end
    end
  end
  %scatter(starts(1,:),starts(2,:),'r');
  %scatter(ends(1,:),ends(2,:),'b');

  return;
end

function [starts, ends, angles] = get_pairs(pt, direct, mindist, maxdist, edge, tangent, convex, precision)

  conv = convex(pt(2),pt(1));
  tang = tangent(pt(2),pt(1));

  candidates = get_line(pt, direct, mindist, maxdist);
  scatter(candidates(1,:),candidates(2,:),'y')
  [candidates, tangs, convs] = get_edges(edge, tangent, convex, candidates)
  scatter(candidates(1,:),candidates(2,:),'k')

  valids = logical(zeros(1,length(convs)));
  for i=1:length(convs)
    t1(i) = check_tangents(tang,tangs(i),precision);
    t2(i) = check_convexity(pt,candidates(:,i),conv,convs(i),precision);
    valids(i) = (check_tangents(tang,tangs(i),precision) && check_convexity(pt,candidates(:,i),conv,convs(i),precision));
  end
  valids

  npairs = sum(valids);
  
  starts = zeros(4,npairs);
  ends = zeros(4,npairs);
  angles = ones(1,npairs);
  if(npairs>0)
    starts(1,:) = pt(1);
    starts(2,:) = pt(2);
    starts(3,:) = tang;
    starts(4,:) =  conv;

    ends(1:2,:) = candidates(:,valids); 
    ends(3,:) = tangs(valids);
    ends(4,:) = convs(valids);

    angles = angles*direct;
  end

  return;
end

function [pts, tangs, convs] = get_edges(edges, direct, convex, pts)

  [h,w] = size(edges);

  pts = pts(:,all(pts>0,1));
  pts = pts(:,(pts(1,:)<=w));
  pts = pts(:,(pts(2,:)<=h));

  indxs = sub2ind([h,w],pts(2,:),pts(1,:));
  isedge = edges(indxs);

  pts = pts(:,isedge);
  indxs = indxs(isedge);

  tangs = direct(indxs);
  convs = convex(indxs);

  return;
end

function pts = get_line(pt1, angl, mindist, maxdist)

  [minx,miny]=pol2cart(angl,mindist);
  [maxx,maxy]=pol2cart(angl,maxdist);

  if(abs(maxx-minx)>abs(maxy-miny))
    xpts = round(minx:sign(maxx-minx):maxx);
    ypts = round(tan(angl)*xpts);
  else
    ypts = round(miny:sign(maxy-miny):maxy);
    xpts = round(ypts/tan(angl));
  end

  pts = [xpts+pt1(1); ypts+pt1(2)];

  [x,y]=pol2cart(angl,maxdist);
  plot(pt1(1)+[0 x],pt1(2)+[0 y],'-or');
  hold on
  plot(xpts+pt1(1),ypts+pt1(2),'og');
  plot([minx maxx]+pt1(1),[miny maxy]+pt1(2),'ob');

  return;
end

function isvalid = check_tangents(angl1, angl2, precision)

  isvalid = (abs(abs(angl2 - angl1) - pi/2) < precision) ;

  return;
end

function isvalid = check_convexity(pt1, pt2, conv1, conv2, precision)

  isvalid = true;

  talph = tan(conv1);
  center = zeros(1,2);

  % position of the point equidistant to pt1 and pt2 (circle center) on the line defined by the curvature of pt1
  % center of the best-fit circle defined by both points
  if(abs(conv1)==pi/2)
    center(1,1) = pt1(1);
    center(1,2) = (pt1(2)^2 - pt2(2)^2 - (pt1(1) - pt2(1))^2) / (2*(pt1(2) - pt2(2)));
  else
    center(1,1) = (pt1(2)^2 - pt1(1)^2 + pt2(1)^2 + pt2(2)^2 + 2*pt2(2)*pt1(1)*talph - 2*pt2(2)*pt1(2) - 2*pt1(1)*pt1(2)*talph)/(2*(pt2(1) - pt1(1) + pt2(2)*talph - pt1(2)*talph));
    center(1,2) = (center(1,1) - pt1(1))*talph + pt1(2);
  end

  alpha1 = atan2(pt1(2)-center(1,2), pt1(1)-center(1,1)) + pi;
  alpha2 = atan2(pt2(2)-center(1,2), pt2(1)-center(1,1)) + pi;

  if(alpha1<-pi)
    alpha1 = alpha1 + 2*pi;
  elseif(alpha1>pi)
    alpha1 = alpha1 - 2*pi;
  end
  if(alpha2<-pi)
    alpha2 = alpha2 + 2*pi;
  elseif(alpha2>pi)
    alpha2 = alpha2 - 2*pi;
  end

  if(conv1<-pi/2 & alpha1>pi/2)
    alpha1 = alpha1 - 2*pi;
  elseif (conv1>pi/2 & alpha1<-pi/2)
    alpha1 = alpha1 + 2*pi;
  end
  if(conv2<-pi/2 & alpha2>pi/2)
    alpha2 = alpha2 - 2*pi;
  elseif (conv2>pi/2 & alpha2<-pi/2)
    alpha2 = alpha2 + 2*pi;
  end

  isvalid = ((abs(conv1-alpha1) < precision) & (abs(conv2-alpha2) < precision));

  conv1
  alpha1
  conv2
  alpha2

  plot([pt1(1) center(1,1)],[pt1(2) center(1,2)],'-^b')
  hold on
  plot([pt2(1) center(1,1)],[pt2(2) center(1,2)],'-or')
  %plot(center(1),center(2),'og')
  %if(pt1(1)<pt2(1))
  %plot([pt1(1):pt2(1)],([pt1(1):pt2(1)]-pt1(1))*talph + pt1(2),'k');
  %else
  %plot([pt1(1):-1:pt2(1)],([pt1(1):-1:pt2(1)]-pt1(1))*talph + pt1(2),'k');
  %end

  theta = [0:2*pi/60:2*pi];
  rads = ones(1,61)*sqrt(sum((pt1-center').^2));
  [x,y]=pol2cart(theta,rads);
  plot(x+center(1),y+center(2),'y');
  [x,y]=pol2cart(conv2,rads(1));
  plot([pt2(1) x+pt2(1)],[pt2(2) y+pt2(2)],'k');

  return;
end

function convex = compute_convexity(edge, direct, nconv)

  [h,w] = size(edge);
  convex = zeros([h w]);
  
  xfilt = [1:nconv] - ceil(nconv/2);

  minbound = -xfilt(1,1);
  maxbound = xfilt(1,end);

  xfilt = repmat(xfilt,nconv,1);
  yfilt = rot90(xfilt);

  symimg = zeros([h w]+minbound+maxbound);
  symimg(minbound+1:end-maxbound,minbound+1:end-maxbound) = edge;

  haxe = find(any(edge,2)');
  waxe = find(any(edge,1));

  for i=haxe
    for j=waxe
      if(edge(i,j))
        cfilt = (yfilt>=(xfilt*tan(direct(i,j))));
        cfilt = cfilt*2 -1;
        convex(i,j) = sum(sum(cfilt.*symimg(i:i+nconv-1,j:j+nconv-1)));
      end
    end
  end

  return
end
