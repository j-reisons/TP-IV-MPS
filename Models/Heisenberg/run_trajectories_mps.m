function run_trajectories_mps(N,U,G,D_max,N_traj)
%%
J = 1;
dt = 0.01;

T_transient = 100; % Time of transient (before averaging)
T_averaging = 3*N; % Time during which averaging occurs
T = T_transient + T_averaging; %Total evolution time
i_cut = round(T_transient/dt);

time = linspace(0,T,floor(T/dt)+1);
[U_odd,U_even] = HeisenbergOpen_U_O2(N,J,U,G,dt);

filename = ['trajectories_MPS','_N',strrep(num2str(N),'.',',') ,'_U',strrep(num2str(U),'.',',')...
    ,'_G',strrep(num2str(G),'.',','),'_Traj',num2str(N_traj),'.mat'];

%rng shuffle
stream = RandStream('mt19937ar','Seed',5489); % MATLAB's start-up settings
RandStream.setGlobalStream(stream);

%%
S_Z =[
    [1 , 0]
    [0 ,-1]
    ];

S_plus = sparse([
    [0,1]
    [0,0]
    ]);

S_minus = S_plus.';
%%
Profiles_all = zeros(N_traj,N);
Currents_all = [];

parfor traj = 1:N_traj
State = random_mps(N,D_max,2);
Profile = zeros(1,N);

for i = 1 : length(time);
        %%%%%%%%%%Record after transient%%%%%%%%
if i > i_cut
    for j = 1:N
        Sz_State = State;
        Sz_State{j} = contract(State{j},3,S_Z,2);
        Profile(j) = Profile(j) + real(braket(Sz_State,State));
    end
end
        %%%%%%%%%%Stochastic evolution%%%%%%%%%%
    S1 = State;
    S1{1} = contract(State{1},3,S_plus,2);
    P1 = real(G*dt*braket(S1,S1));
    
    S2 = State;
    S2{N} = contract(State{N},3,S_minus,2);
    P2 = real(G*dt*braket(S2,S2));

    e = rand();
    
    if e <= P1
        State = sweep(S1,1);
    elseif e <= (P2 + P1) && e >= P1
        State = sweep(S2,1);
    else
        State = apply(U_odd,State);
        State = sweep(State,1);
        State = sweep(State,-1,D_max);
        State = apply(U_even,State);
        State = sweep(State,-1);
        State = sweep(State,1,D_max);
        State = apply(U_odd,State);
        State = sweep(State,1);
        State = sweep(State,-1,D_max);
    end
end
Profile = Profile/(length(time)- i_cut);
Profiles_all(traj,:) = Profile;
Currents_all(traj) = (Profile(N) + 1)/2;
end

Mean_Profile = mean(Profiles_all,1);
Mean_Current = mean(Currents_all);

save(filename,'Profiles_all','Currents_all','Mean_Profile','Mean_Current');
end