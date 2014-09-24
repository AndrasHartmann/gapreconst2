function [stable] = checkstable(A);
% function [stable] = checkstable(A);
% simple stability check Using Durbin recursion as in [1]
%
% Return:
%   stable: true if the polynom A holds coeffitients of a stable IIR filter
%
%References
% [1] https://ccrma.stanford.edu/~jos/fp/Testing_Filter_Stability_Matlab.html

N = length(A)-1; % Order of A(z)
stable = 1;      % stable unless shown otherwise
A = A(:);        % make sure it's a column vector
for i=N:-1:1
  rci=A(i+1);
  if abs(rci) >= 1
    stable=0;
    return;
  end
  A = (A(1:i) - rci * A(i+1:-1:2))/(1-rci^2);
  % disp(sprintf('A[%d]=',i)); A(1:i)'
end
