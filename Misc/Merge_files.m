function Merge_files(merged_filenames_cell)

Merged_Currents_all = [];
Merged_Cut_data_all = [];
Merged_D_of_t_all = [];
Merged_Profiles_all = [];

%%%
%%%
for i = 1:length(merged_filenames_cell)
    load(merged_filenames_cell{i})
    Merged_Currents_all = [Merged_Currents_all;Currents_all];
    Merged_Cut_data_all = [Merged_Cut_data_all;Cut_data_all];
    Merged_D_of_t_all = [Merged_D_of_t_all;D_of_t_all];
    Merged_Profiles_all = [Merged_Profiles_all;Profiles_all];
end

Profiles_all = Merged_Profiles_all;
Currents_all = Merged_Currents_all;
Cut_data_all = Merged_Cut_data_all;
D_of_t_all = Merged_D_of_t_all;

Mean_Profile = mean(Profiles_all,1);
Mean_Current = mean(Currents_all);

N_traj = length(Currents_all);

filename = ['trajectories_MPS','_N',strrep(num2str(N),'.',',') ,'_U',strrep(num2str(U),'.',',')...
    ,'_G',strrep(num2str(G),'.',','),'_tol',num2str(tolerance),'_D',num2str(D_limit)...
    ,'_Traj',num2str(N_traj),'.mat'];

save(filename,'N','U','G','tolerance','D_limit','Profiles_all'...
    ,'Currents_all','Mean_Profile','Mean_Current','Cut_data_all','D_of_t_all');

end