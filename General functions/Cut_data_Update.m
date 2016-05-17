function Cut_data = Cut_data_Update(error_v,D_v,Cut_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Cut_data(1) = Cut_data(1) + sum(error_v);
Cut_data(2) = max([Cut_data(2),error_v]);
Cut_data(3) = max([Cut_data(3),D_v]);
end

