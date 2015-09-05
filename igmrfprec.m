function [Q,W,q]=igmrfprec(sz,order)
% IGMRFPREC Constructs a precision matrix for a 1:st or 2:nd order IGMRF
%
%   Q = igmrfprec(sz,order)
%   [Q,W] = igmrfprec(sz,order)
%   [Q,W,q] = igmrfprec(sz,order)
%
%   sz : The grid/image size,  sz = [m,n] for an m-by-n image
%   order : The order of the increments, 1 or 2.
%   Q  : A sparse precision matrix
%   W  : The increment matrix: Q=W'*W
%   q  : The precision coefficients in the format used
%        for input to gmrfprec.
%        The Q from gmrfprec(sz,q) differs from the Q from igmrfprec
%        only for pixels near the image boundary.
%
% For order==1, the increments are:
%   [-1 1] and [-1; 1]
% The resulting field density is invariant to addition of a constant.
%
% For order==2, the increment for all interior grid points is
%   [0 1 0; 1 -4  1; 0 1 0]
% The increments for the horizontal and vertical edges are
%   [1 -2 1] and [1; -2; 1]
% The increments for the four corners are
%   [-1 1; 1 -1]
% The boundary increments are motivated by assuming that the second
% derivative perpendicular to the boundary is zero.
% The resulting field density is invariant to addition of an arbitrary
% plane.
%
% See also: gmrfprec

% $Id: igmrfprec.m 2201 2005-11-14 01:14:25Z finn $

n = prod(sz);

if (order==1)
  N = prod(sz-[1,0])+prod(sz-[0,1]);
  
  % Vertical increments
  [I,J] = ndgrid(1:sz(1)-1,1:sz(2));
  I = I(:); J = J(:);
  M1 = length(I);
  II = [1:M1,1:M1];
  JJ = [(I  +sz(1)*(J-1));...
        (I+1+sz(1)*(J-1))];
  KK = [ones(M1,1);-ones(M1,1)];
  M = M1;
  
  % Horizontal increments
  [I,J] = ndgrid(1:sz(1),1:sz(2)-1);
  I = I(:); J = J(:);
  M2 = length(I);
  II = [II,M+(1:M2),M+(1:M2)];
  JJ = [JJ;...
        (I+sz(1)*(J  -1));...
        (I+sz(1)*(J+1-1))];
  KK = [KK;ones(M2,1);-ones(M2,1)];
  M = M+M2;

  W = sparse(II,JJ,KK,N,n);
  Q = W'*W;

  if (nargout>=3)
    Q_mini = igmrfprec([3,3],1);
    q = full(icolstack(Q_mini(:,5),[3,3]));
  end
elseif (order==2)
  % N = prod(sz-1)+sum(sz-2)+4;

  if (sz(1)>=3) & (sz(2)>=3)
    % Internal increments
    [I,J] = ndgrid(2:sz(1)-1,2:sz(2)-1);
    I = I(:); J = J(:);
    M = length(I);
    II = kron(ones(5,1),(I+sz(1)*(J-1)));
    JJ = [(I-1+sz(1)*(J  -1));...
          (I+1+sz(1)*(J  -1));...
          (I  +sz(1)*(J-1-1));...
          (I  +sz(1)*(J+1-1));...
          (I  +sz(1)*(J  -1))];
    KK = [ones(M*4,1);-4*ones(M,1)];
    M_ = M;
  else
    II = [];
    JJ = [];
    KK = [];
    M_ = 0;
  end
  
  if (sz(1)>=3)
    % Vertical edge increments
    [I,J] = ndgrid(2:sz(1)-1,[1,sz(2)]);
    I = I(:); J = J(:);
    M = length(I);
    II = [II;kron(ones(3,1),(I+sz(1)*(J-1)))];
    JJ = [JJ;...
          (I-1+sz(1)*(J-1));...
          (I  +sz(1)*(J-1));...
          (I+1+sz(1)*(J-1))];
    KK = [KK;ones(M,1);-2*ones(M,1);ones(M,1)];
    M_ = M_+M;
  end
  
  if (sz(2)>=3)
    % Horizontal edge increments
    [I,J] = ndgrid([1,sz(1)],2:sz(2)-1);
    I = I(:); J = J(:);
    M = length(I);
    II = [II;kron(ones(3,1),(I+sz(1)*(J-1)))];
    JJ = [JJ;...
          (I+sz(1)*(J-1-1));...
          (I+sz(1)*(J  -1));...
          (I+sz(1)*(J+1-1))];
    KK = [KK;ones(M,1);-2*ones(M,1);ones(M,1)];
    M_ = M_+M;
  end

  if (sz(1)>=2) & (sz(2)>=2)
    % Corner increments
    [I,J] = ndgrid(unique([1,sz(1)-1]),unique([1,sz(2)-1]));
    I = I(:); J = J(:);
    M = length(I);
    II = [II;kron(ones(4,1),([1;sz(1);1;sz(1)]+...
                             sz(1)*([1;1;sz(2);sz(2)]-1)))];
    JJ = [JJ;...
          (I  +sz(1)*(J  -1));...
          (I+1+sz(1)*(J+1-1));...
          (I+1+sz(1)*(J  -1));...
          (I  +sz(1)*(J+1-1))];
    KK = [KK;[ones(M*2,1);-ones(M*2,1)]];
    M_ = M_+M;
  end

  W = sparse(II,JJ,KK,max(II),n);
  Q = W'*W;
  
  if (nargout>=3)
    Q_mini = igmrfprec([5,5],2);
    q = full(icolstack(Q_mini(:,13),[5,5]));
  end
else
  error(sprintf('Orders >2 not implemented.',order))
end


% -------------------
function x=icolstack(y,sz);
% ICOLSTACK Invert column stacking of an image 
%
%  x=icolstack(y,sz)
%
%   sz: [m,n], size(x) = [sz,size(y,2)]
%       or
%       sz = size(x)
%
% SEE ALSO colstack

% $Id: icolstack.m 3318 2007-04-04 20:31:03Z finn $

x=reshape(full(y),[sz(1),sz(2),size(y,2)]);
