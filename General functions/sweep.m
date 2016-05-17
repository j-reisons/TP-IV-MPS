function [mps_out,error_vector,D_vector] = sweep(mps_in,direction,varargin)
% varargin is either a precision or a maximal bond dimension
% For precision / tolerance cut, a maximal bond dimension of 500 is
% enforced

mps_out = mps_in;
N = length(mps_in);
error_vector = zeros(1,N);
D_vector = zeros(1,N);
switch direction
    case 1
        for i = 1:N
            [mps_out,error,D] = L_can(mps_out,i,varargin{:});
            error_vector(i) = error;
            D_vector(i) = D;
        end
        
    case -1
        for i = N:-1:1
            [mps_out,error,D] = R_can(mps_out,i,varargin{:});
            error_vector(i) = error;
            D_vector(i) = D;
        end
        
end