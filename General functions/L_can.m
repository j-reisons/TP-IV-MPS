function [Mleft,error,D_max] = L_can(mps,site,varargin)
% Left canonizes (and compresses) site of mps and throws the SV' to the right
% varargin is empty, or tolerance and D_limit
% If empty, canonization without compression
% If tolerance and D_lmit is provided compresses to tolerance or D_limit
% Last site normalization is thrown out
a = 0;
if ~isempty(varargin)
    tolerance = varargin{1};
    D_limit = varargin{2};
    a = 1;
end

N = length(mps);
Mleft = mps;

% Middle sites    
work = Mleft{site};
s_w = size(work);
work = permute(work,[1 3 2]);
work = reshape(work,[s_w(3)*s_w(1),s_w(2)]);
try
    [U,S,V] = svd(work,'econ');
catch
    work = work + rand(size(work))*1E-10;
    [U,S,V] = svd(work,'econ');
end
s_s = size(S,1);

switch a
    case 0
    S = S /sqrt(trace(S*S'));
    D_max = size(S,1);
    error = 0;
    case 1
    S_2 = S*S';
    S = S /sqrt(trace(S_2));
    S_2 = diag(S*S');
    cut = 0;
    sum = 0;
    while 1 - sum > tolerance && cut < D_limit && cut < length(S_2)
        cut = cut +1;
        sum = sum + S_2(cut);
    end
    error = 1 - sum;
    D_max = cut;
    U = U(:,1:D_max);
    S = S(1:D_max,1:D_max);
    S = S /sqrt(trace(S*S'));
    V = V(:,1:D_max);
end
s_u = size(U);
A = reshape(U,[s_w(1) s_w(3) s_u(2)]);
A = permute(A,[1 3 2]);
Mleft{site} = A;
SV = S*V';
if site ~= N
    Mleft{site+1} = contract(SV,2,Mleft{site+1},1);
else
    Mleft{site} = Mleft{site}*sign(SV); % Fixing +- 1 phase
end

end