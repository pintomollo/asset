function [clusts] = stsc(D, CLUSTER_NUM_CHOICES, neighbor_num)
% Self-Tuning Spectracl Clustering
%
% c.f: http://www.vision.caltech.edu/lihi/Demos/SelfTuningClustering.html
%
%%  Code by Lihi Zelnik-Manor (2005)

  [D_LS,A_LS] = scale_dist(D, floor(neighbor_num/2)); %% Locally scaled affinity matrix

  [clusts_RLS, rlsBestGroupIndex, qualityRLS] = cluster_rotate(A_LS,CLUSTER_NUM_CHOICES,0,1)

  %keyboard

  tmp_clusts = clusts_RLS{rlsBestGroupIndex};
  clusts = zeros(size(D, 1), 1);

  for i = 1:numel(tmp_clusts)
    clusts(tmp_clusts{i}) = i;
  end

  return;
end

function [clusts,best_group_index,Quality,Vr] = cluster_rotate(A,group_num,fig,method)

%% cluster by rotating eigenvectors to align with the canonical coordinate
%% system
%%
%%   [clusts,best_group_index,Quality,Vr] = cluster_rotate(A,group_num,method,fig)
%%  
%%  Input:
%%        A = Affinity matrix
%%        group_num - an array of group numbers to test
%%                    it is assumed to be a continuous set
%%        fig - Figure to display progress. set to 0 if no display is
%%              desired
%%        method - 1   gradient descent 
%%                 2   approximate gradient descent
%%        
%%  Output:
%%        clusts - a cell array of the results for each group number
%%        best_group_index - the group number index with best alignment quality
%%        Quality = the final quality for all tested group numbers
%%        Vr = the rotated eigenvectors
%%
%%
%%  Code by Lihi Zelnik-Manor (2005)
%%
%%


  if( nargin < 2 )
      group_num = [2:6];
  end
  if( nargin < 3 )
      fig = 0;
  end
  if( nargin < 4 )
      method = 1;  %% method to calculate cost gradient. 1 means true derivative
                   %% change to any other value to estimate fradient numerically
  end
  group_num = sort(group_num);
  group_num = setdiff(group_num,1);

  %%% obtain eigenvectors of laplacian of affinity matrix
  %tic; 
  nClusts = max(group_num);
  [V,evals] = evecs(A,nClusts); 
  %ttt = toc;
  %disp(['evecs took ' num2str(ttt) ' seconds']);

  %%%%%% Rotate eigenvectors
  %clear clusts;
  Vcurr = V(:,1:group_num(1));
  for g=1:length(group_num),
      %%% make it incremental (used already aligned vectors)
      if( g > 1 )
          Vcurr = [Vr{g-1},V(:,group_num(g))];
      end
      [clusts{g},Quality(g),Vr{g}] = evrot(Vcurr,method);
  end
  i = find(max(Quality)-Quality <= 0.001);
  best_group_index = i(end);
end

function [V,ss,L] = evecs(A,nEvecs)

%% calculate eigenvectors, eigenvalues of the laplaican of A
%%
%%   [V,ss,L] = evecs(A,nEvecs)
%%  
%%  Input:
%%        A = Affinity matrix
%%        nEvecs = number of eigenvectors to compute
%%        
%%  Output:       
%%        V = eigenvectors
%%        ss = eigenvalues
%%        L = Laplacian
%%
%%
%%  Code by Lihi Zelnik-Manor (2005)
%%
%%



  %%%%%%%% Compute the Laplacian
  %tic;
  npix = size(A,1);
  useSparse = issparse(A);
  dd = 1./(sum(A)+eps);
  dd = sqrt(dd);
  if(useSparse)
      DD = sparse(1:npix,1:npix,dd);
  else
      DD = diag(dd);
  end
  L = DD*A*DD;
  %ttt = toc;
  % disp(['Laplacian computation took ' num2str(ttt) ' seconds']);


  %%%%%%% Compute eigenvectors
  %tic;
  if (useSparse)
      opts.issym = 1;
      opts.isreal = 1;
      opts.disp = 0;
      [V,ss] = eigs(L,nEvecs,1,opts);
  %     [VV,ss]=svds(L,nClusts,1,opts);
  else
      [V,ss] = svd(L);
      V = V(:,1:nEvecs);    
  end
  ss = diag(ss);
  ss = ss(1:nEvecs);
  %ttt = toc;
  % disp(['eigenvectors computation took ' num2str(ttt) ' seconds']);
end
