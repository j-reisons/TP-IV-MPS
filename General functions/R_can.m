function Mright = R_can(mps,site,varargin)
% Right canonizes (and compresses) site of mps and throws the US to the left.
% If D_max is provided, also compresses the bond to D_max bond dimension
% First site normalization is thrown out

a = 0; % 0 is no cut, 1 is D_max cut, 2 is tolerance cut

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
Mright = mps;

work = Mright{site};
s_w = size(work);
work = reshape(work,[s_w(1),s_w(2)*s_w(3)]);
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
    case 1
        if size(S,1) > D_max
    U = U(:,1:D_max);
    S = S(1:D_max,1:D_max);
    S = S /sqrt(trace(S*S'));
    V = V(:,1:D_max);
        end
    case 2
    S_2 = S*S';
    S = S /sqrt(trace(S_2));
    S_2 = diag(S_2);
    cut = 1;
    sum = 0;
    while 1 - sum > tolerance
        sum = sum + S_2(cut);
        cut = cut +1;
    end
    D_max = cut;
    U = U(:,1:D_max);
    S = S(1:D_max,1:D_max);
    S = S /sqrt(trace(S*S'));
    V = V(:,1:D_max);
end
s_v = size(V');
B = reshape(V',[s_v(1),s_w(2),s_w(3)]);
Mright{site} = B;
US = U*S;
if site ~=1
    Mright{site-1} = contract(Mright{site-1},2,US,1);
    Mright{site-1} = permute(Mright{site-1},[1,3,2]);
else
    Mright{site} = Mright{site}*sign(US); % Fixing +- 1 phase
end

end

