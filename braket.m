function sprod = braket(mps_1,mps_2)
% Returns the scalar product < MPS_1 | MPS_2 > that runs in polynomial time.
% The inputs are cell arrays corresponding to MPS decompositions,
% notice that the bond dimensions of the MPS does not have to be the same.

N = length(mps_1);
D_S = size(mps_1{1},3);
sprod = 0;

for sigma = 1:D_S
    sprod = sprod + mps_1{1}(:,:,sigma)'*mps_2{1}(:,:,sigma);
end

for i = 2:N
    sum = 0;
    for sigma= 1:D_S
        sum = sum + mps_1{i}(:,:,sigma)'*sprod*mps_2{i}(:,:,sigma);
    end
    sprod = sum;
end
end