function [M,x] = grid_composite(varargin)
% Creates the differentiation matrix for a composite basis consisting of n basis segments
% - Enter n pairs of arguments;
%       x_i - grid, column vector
%       M_i - differentiation matrix corresponding to x_i
% - bm - method for joining at boundaries; identified is number of inputs is odd
%           0 - value of derivative at end points is average of either side (default)
%           1 - value of derivative at endpoints is continuous, set by left value
%           2 - value of derivative at endpoints is continuous, set by right value
%
% ----------------------------------------------------------------------------
% Note: You may join an arbitrary number of basis segments though some
%       choices make more sense than others, e.g. including a Laguerre
%       basis segment in the middle of a composite basis may lead to
%       incorrect results as this basis assumes the domain is semi-infinite
%       in one direction.
%
% Note: The endpoints of each segment much correspond to those of the
%       surrounding segments, e.g. the segments [1 2 3] and [3 4 5] may be
%       joined as 3 is in both. Conversely, a type 1 Chebyshev basis (see
%       diff_mat.m) on [-1 0] cannot be joined to a type 2 Chebyshev basis 
%       on [0 1] as the first segment does not include 0.
% ----------------------------------------------------------------------------

if nargin/2==floor(nargin/2)
    n = nargin/2;
    bm = 0;
else
    n = (nargin-1)/2;
    bm = varargin{nargin};
end

% define x, check endpoints

x = varargin{1};
N = length(x)*ones(1,n);

for in = 2:n
	if varargin{2*in-1}(1) ~= x(end); error('Endpoint mismatch'); end
	x = [x; varargin{2*in-1}(2:end)];
    N(in) = length(varargin{2*in-1});
end

M = zeros(sum(N)+1-n);

for in = 1:n
    
    i1 = sum(N(1:in-1))+2-in;
    i2 = sum(N(1:in))-in+1;
    
    if bm == 0
        if in > 1
            M(i1,:) = M(i1,:)/2; M(i1,i1:i2) = M(i1,i1:i2)+varargin{2*in}(1,:)/2;
            M(i1+1:i2,i1:i2) = varargin{2*in}(2:end,:);
        else
            M(i1:i2,i1:i2) = varargin{2*in};
        end
    end
    
    if bm == 1
        if in > 1
            M(i1+1:i2,i1:i2) = varargin{2*in}(2:end,:);
        else
            M(i1:i2,i1:i2) = varargin{2*in};
        end
    end
    
	if bm == 2
        if in > 1
            M(i1,:) = 0;
            M(i1:i2,i1:i2) = varargin{2*in};
        else
            M(i1:i2,i1:i2) = varargin{2*in};
        end
	end
    
end

end

