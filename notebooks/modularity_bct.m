%%

cd '/home/finc/Dropbox/Projects/LearningBrain/data/neuroimaging/'
%%

M = load('correlation_matrices_static_lp_0_08.mat');
dual = M.correlation_matrices_all;


%%
MOD = zeros(size(dual,3), 4);
nrep = 200;
n = 1;
dual(isnan(dual(:))) = 0;

%%


for sub = 1: size(dual,1)
    sub
  
    for ses = 1: size(dual,2)
        ses
        for con = 1: size(dual,3)
            qt = 0;
            for rep = 1: nrep
                [Mod, Q] = community_louvain(squeeze(dual(sub,ses,con,:,:)),[],[],'negative_asym');
                    if Q > qt
                    qt = Q;
                    end
            end
                    MOD(n, 1) = sub; % subjects
                    MOD(n, 2) = ses; % sessions
                    MOD(n, 3) = con; % sessions
                    MOD(n, 4) = qt; % sessions
                    n = n + 1;
end
end
end
MOD