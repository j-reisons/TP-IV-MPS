function [Mleft,error,D_max] = L_can(mps,site,varargin)
% Left canonizes (and compresses) site of mps and throws the SV' to the right
% varargin is empty, D_max, or tolerance
% If empty, canonization without compression
% If D_max is provided, compresses the bond to D_max bond dimension
% If tolerance is provided, compresses to tolerance or D_max of 500
% Last site normalization is thrown out
a = 0;
if ~isempty(varargin)
    if mod(varargin{1},1) == 0
        D_max = varargin{1};
        a = 1;
    elseif varargin{1} < 1 && varargin{1} > 0
        tolerance = varargin{1};
        a = 2;
    end
end

N = length(mps);
Mleft = mps;

% Middle sites    
work = Mleft{site};
s_w = size(work);
work = permute(work,[1 3 2]);
work = reshape(work,[s_w(3)*s_w(1),s_w(2)]);
s_w2 = size(work);
try
    [U,S,V] = svd(work,'econ');
catch
    work = work + rand(s_w2)*1E-12;
    [U,S,V] = svd(work,'econ');
end
s_s = size(S,1);

switch a
    case 0
    S = S /sqrt(trace(S*S'));
    D_max = size(S,1);
    error = 0;
    case 1
        if size(S,1) > D_max
            S = S /sqrt(trace(S*S'));
            U = U(:,1:D_max);
            S = S(1:D_max,1:D_max);
            error = 1 - trace(S*S');
            S = S /sqrt(trace(S*S'));
            V = V(:,1:D_max);
        else
            error = 0;
        end
    case 2
    S_2 = S*S';
    S = S /sqrt(trace(S_2));
    S_2 = diag(S*S');
    cut = 0;
    sum = 0;
    while 1 - sum > tolerance && cut < 500 && cut < length(S_2)
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