function n = sparsity(x,threshold)
% n = sparsity(x,threshold) returns the number of elements in x that
% contain most of the weight of the vector.
% $Id: sparsity.m 326 2008-06-25 18:34:35Z mpf $
if nargin < 2,
  threshold = 0.9995;
else
  threshold = max(0,min(threshold,1));
end

x = sort(abs(x),'descend');
if sum(x) == 0
  n = 0;
else
  x = x / sum(x);
  x = cumsum(x);
  n = min(find(x >= threshold));
end
