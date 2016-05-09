function run_trajectories(N,U,G,N_traj)

%rng shuffle
stream = RandStream('mt19937ar','Seed',5489); % MATLAB's start-up settings
RandStream.setGlobalStream(stream);

%% Simulation parameters

% H = Sum over i = [0:1:N-1] { J*( S_X(i)*S_X(i+1) + S_Y(i)*S_Y(i+1) )
% + U( S_Z(i)*S_Z(i+1) ) }

filename = ['N',strrep(num2str(N),'.',',') ,'_U',strrep(num2str(U),'.',',')...
    ,'_G',strrep(num2str(G),'.',','),'_trajectories.mat'];

J = 1; 

dt = 10^-3; % Time increment

T = 2*(1/J)*N; %Total evolution time
%% Pauli and co.

S_X = sparse([
    [0 , 1]
    [1 , 0]
    ]);

S_Y = sparse([
    [0 , -1i]
    [1i , 0]
    ]);

S_Z = sparse([
    [1 , 0]
    [0 ,-1]
    ]);

S_plus = sparse([
    [0,1]
    [0,0]
    ]);

S_minus = S_plus.';

I = speye(2);


%% Building the Hamiltonian, S_plus_1 , S_minus_N , S_Z_indexed matrices

H = sparse(0);

for i = 1:N-1
    O_X = 1;
    O_Y = 1;
    O_Z = 1;
    for j = 1:N
        if j == i || j == i+1
            O_X = kron(O_X,S_X);
            O_Y = kron(O_Y,S_Y);
            O_Z = kron(O_Z,S_Z);
        else
            O_X = kron(O_X,I);
            O_Y = kron(O_Y,I);
            O_Z = kron(O_Z,I);
        end
    end
    H = H + J*O_X + J*O_Y + U*O_Z ;
end

S_plus_1 = kron(S_plus,speye(2^(N-1)));

S_minus_N = kron(speye(2^(N-1)),S_minus);

S_Z_indexed = cell(N);

for i=1:N
    S_Z_i = kron(speye(2^(i-1)),S_Z);
    S_Z_i = kron(S_Z_i,speye(2^(N-i)));
    S_Z_indexed{i} = S_Z_i;
end

%% Effective hamiltonian, M(x,t)

H_eff = H -(((1i)*G)/2)*(S_plus_1'*S_plus_1) - (((1i)*G)/2)*(S_minus_N'*S_minus_N);

M = zeros(N,round(T/dt) + 1);

%%
parfor k = 1:N_traj
    
    State = rand(2^N,1);
    State = State / norm(State);
    
    M_partial = zeros(N,round(T/dt) + 1);
    
    for time = 0:dt:T
        S1 = sqrt(G*dt)*(S_plus_1*State);
        P1 = norm(S1)^2;
        
        S2 = sqrt(G*dt)*(S_minus_N*State);
        P2 = norm(S2)^2;
        
        e = rand;
        
        %Doing the stochastic evolution
        
        if e <= P1
            State = S1/sqrt(P1);
        elseif e <= (P2 + P1) && e >= P1
            State = S2/sqrt(P2);
        else
            State = State - (1i)*dt*(H_eff*State);
            State = State/norm(State);
        end
        
        % Computing M(x,t)
        step = round(time/dt) + 1;
        for x=1:N
            M_partial(x,step) = State'*((S_Z_indexed{x})*State);
        end
        
    end
    M = M + M_partial
end

M = M/N_traj;


a = round(length(M)/2);
b = length(M);
M_cut = M(:,a:b);
M_mean = mean(M_cut,2);

J_N = (M_mean(N) + 1)/2;

save(filename,'M_mean','J_N');