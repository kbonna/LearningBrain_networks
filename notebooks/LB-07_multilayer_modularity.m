%W = squeeze(readNPY('/home/finc/Dropbox/Projects/LearningBrain/data/neuroimaging/LearningBrain_dFC_matrices.npy'));

large = load('/home/finc/Dropbox/Projects/LearningBrain/data/neuroimaging/LearningBrain_dFC_matrices_dualnback.mat');


%%
W = squeeze(large.correlation_matrices_dFC);
W = W(1:10, :, :, :, :);

%%
%WT = zeros(46, 4, 20, 264, 264);
WT = zeros(10, 4, 20, 264, 264);


%%
for i = 1: length(W(:,1,1,1,1))
    for j = 1: length(W(1,:,1,1,1)) 
        for k = 1: length(W(1,1,:,1,1))
            T = threshold_absolute(squeeze(W(i,j,k,:,:)), 0);
            WT(i,j,k,:,:) = T;
            %WT(i,j,k,:,:) = squeeze(W(i,j,k,:,:).* squeeze(W(i,j,k,:,:)) > 0);
            
        end
    end
end


%% calculate modularity for one session

% ------ create cell array filled by matrices

A = cell(1,20);

for k = 1: length(WT(1,1,:,1,1))
    A{k} = squeeze(WT(1,1,k,:,:));
end

%% create multilayer network

omega = 1
gamma = 1
N = length(A{1});                       % size of matrix
T = length(A);                          % no of layers


B = spalloc(N*T, N*T, N*N*T + 2*N*T);   % empty multilayer matrix
twomu = 0;


for s = 1:T 
    k = sum(A{s});
    twom = sum(k);
    twomu = twomu+twom;
    indx = [1:N] + (s-1)*N;
    B(indx,indx) = A{s}-gamma*k'*k/twom;
end

twomu = twomu + 2*omega*N*(T-1);

% --------- create categorical multilayer network
%for d = 1:T
%    B = B + omega*spdiags(ones(N*T,2),[-d*N,N],N*T,N*T);
%end

B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);



%% --- calculate multilayer modules -------------------------------------

        

%%

n_sub = 1

modularity = zeros(n_sub, 4);
flex = zeros(n_sub, 4, 264);
flex_mean = zeros(n_sub, 4);     
allegiance = zeros(n_sub, 4, 264, 20);

%%

gamma = 1;
omega = 1;
rep = 10;

for sub = 1:10
    for ses = 1:4
        A = cell(1,20);

        for k = 1:20
            A{k} = squeeze(WT(sub,ses,k,:,:));
            
            N = length(A{1});                       % size of matrix
            T = length(A);                          % no of layers


            B = spalloc(N*T, N*T, N*N*T + 2*N*T);   % empty multilayer matrix
            twomu = 0;
        end

           
        for s = 1:T 
            n = sum(A{s});
            twom = sum(n);
            twomu = twomu+twom;
            indx = [1:N] + (s-1)*N;
            B(indx,indx) = A{s}-gamma*n'*n/twom;
        end

        twomu = twomu + 2*omega*N*(T-1);
        B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
                                 
         Qb = 0;
         for i = 1 : rep
            clc;
                %fprintf('Subject = %i\n',sub);
           [St,Qt] = genlouvain(B);
           Qt = Qt / twomu;
           if Qt > Qb 
              Qb = Qt;
              Sb = reshape(St, N, T);
           end
         end
         
         f = flexibility(Sb', 'temp');
         
         modularity(sub, ses) = Qb;
         flex(sub, ses, :) = f;
         flex_mean(sub, ses) = mean(f);
         allegiance(sub, ses, :, :) = Sb;
         
         
   end
end





