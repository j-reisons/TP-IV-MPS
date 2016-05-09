function mps_out = sweep(mps_in,direction,varargin)
% varargin is either a precision or a maximal bond dimension

mps_out = mps_in;
N = length(mps_in);
switch direction
    case 1
        for i = 1:N
            mps_out = L_can(mps_out,i,varargin{:});
        end
        
    case -1
        for i = N:-1:1
            mps_out = R_can(mps_out,i,varargin{:});
        end
        
end