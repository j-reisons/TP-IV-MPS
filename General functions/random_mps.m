function [mps] = random_mps(N,D,d,varargin)
%RANDOM_MPS random mps of spin chain length N, bond length D, local H space
%dimension d. Returns Right or left normalized if specified, by default
%left



mps = cell(1,N);
mps{1} = rand(1,D,d);
for i = 2:N-1
   mps{i} = rand(D,D,d);
end
mps{N} = rand(D,1,d);

if ~isempty(varargin)
    direction = varargin{1};
switch direction
    case 1
        mps = sweep(mps,1);
        return
    case -1
        mps = sweep(mps,-1);
        return
end
end
mps = sweep(mps,1);
end

