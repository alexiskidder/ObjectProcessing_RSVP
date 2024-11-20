%% bayes for difference in rsa model fits between binarized and intact stimuli
% AK - 10/21/24
% used output from new_linearmodel.m



%% set up

clear all 

% directories
datadir = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/rdm/" ;% where data is from
modeldir = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/modelRDMS/" ; % where model RDMs are located
old_outdir = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/modelrdm/" ; %where output will be placed
outdir = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/modelrdm/new/" ;

% load 
    % all dimensions 

origbeta_file_all = load(outdir + "zscore_all_participants_orig_all.mat"); % all participants with each dimension for intact stimuli
temp_orig_all = origbeta_file_all.orig_dimall_betas;

binbeta_file_all = load(outdir + "zscore_all_participants_bin_all.mat"); % all participants with each dimension for binarized stimuli
temp_bin_all = binbeta_file_all.bin_dimall_betas;

    % aspcat
origbeta_file_aspcat = load(outdir + "zscore_all_participants_orig_aspcat.mat"); % all participants with asp and cat dimensions for intact stimuli
temp_orig_aspcat = origbeta_file_aspcat.orig_aspcatall_betas;

binbeta_file_aspcat = load(outdir + "zscore_all_participants_bin_aspcat.mat"); % all participants with asp and cat dimensions for binarized stimuli
temp_bin_aspcat = binbeta_file_aspcat.bin_aspcatall_betas;

% where to put betas for all participants

all_asp_orig = [];
all_cat_orig = [];
all_anim_orig = [];
all_asp_bin = [];
all_cat_bin = [];
all_anim_bin = [];

aspcat_asp_orig = [];
aspcat_cat_orig = []; 
aspcat_asp_bin = [];
aspcat_cat_bin = []; 

% put data in correct matrices (order is aspect ratio, animacy, category
% for all models) (order is aspect ratio, category for aspcat)

for part = 1:size(temp_bin_aspcat,2)

    % get participant-specific betas
    part_all_orig = temp_orig_all{part};
    part_all_bin = temp_bin_all{part};
    part_aspcat_orig = temp_orig_aspcat{part};
    part_aspcat_bin = temp_bin_aspcat{part};

    %sort into different dimensions and model types 
    all_asp_orig(part,:) = part_all_orig(1,:);
    all_cat_orig(part,:) = part_all_orig(3,:);
    all_anim_orig(part,:) = part_all_orig(2,:);
    all_asp_bin(part,:) =  part_all_bin(1,:);
    all_cat_bin(part,:) = part_all_bin(3,:);
    all_anim_bin(part,:) = part_all_bin(2,:);

    aspcat_asp_orig(part,:) = part_aspcat_orig(1,:);
    aspcat_cat_orig(part,:) = part_aspcat_orig(2,:);
    aspcat_asp_bin(part,:) =  part_aspcat_bin(1,:);
    aspcat_cat_bin(part,:) = part_aspcat_bin(2,:);

end


%% calculate differences 

% set matrices
all_asp_diff = [];
all_anim_diff = [];
all_cat_diff = [];

aspcat_asp_diff = [];
aspcat_cat_diff = [];

% calculate differences for each participant (orig - bin)

for p = 1:size(all_asp_orig,1)
    for s = 1:size(all_asp_orig,2)
        all_asp_diff(p,s) =  all_asp_orig(p,s) - all_asp_bin(p,s);
        all_anim_diff(p,s) = all_anim_orig(p,s) - all_anim_bin(p,s);
        all_cat_diff(p,s) = all_cat_orig(p,s) - all_cat_bin(p,s);

        aspcat_asp_diff(p,s) = aspcat_asp_orig(p,s) - aspcat_asp_bin(p,s);
        aspcat_cat_diff(p,s) = aspcat_cat_orig(p,s) - aspcat_cat_bin(p,s);
    end
end


%% bayes!

% set parameters
mval = ['mu=0, rscale="medium",nullInterval=c(-0.5,0.5)'];

% set details
for a = 1:5
    if a==1
        file = all_asp_diff;
        outfile = [outdir + 'all_bf_aspdiff_p=20.mat']; 
        disp('starting bayes for all-asp differences')
    elseif a==2
        file = all_anim_diff;
        outfile = [outdir + 'all_bf_animdiff_p=20.mat']; 
        disp('starting bayes for all-anim differences')  
    elseif a == 3 
        file = all_cat_diff;
        outfile = [outdir + 'all_bf_catdiff_p=20.mat']; 
        disp('starting bayes for allasp ratio differences')
    elseif a == 4
        file = aspcat_asp_diff;
        outfile = [outdir + 'aspcat_bf_aspdiff_p=20.mat'];
        disp('starting bayes for aspcat-asp differences')
    elseif a == 5
      file = aspcat_cat_diff;
      outfile = [outdir + 'aspcat_bf_catdiff_p=20.mat'];
      disp('starting bayes for aspcat-cat differences')  
    end

    % stats
    bf = [];
    X = file';
    bf(1,:) = bayesfactor_R_wrapper(X,'args',mval,'returnindex',2);
    
    save(outfile,"bf")

end


%% make files for plotting in python

% get average difference score
mean_all_asp = mean(all_asp_diff);
mean_all_anim = mean(all_anim_diff);
mean_all_cat = mean(all_cat_diff);

mean_aspcat_asp = mean(aspcat_asp_diff);
mean_aspcat_cat = mean(aspcat_cat_diff);

% get bayes info
all_asp_bayes = load(outdir + 'all_bf_aspdiff_p=20.mat');
all_asp_bayes = all_asp_bayes.bf;

all_anim_bayes = load(outdir + 'all_bf_animdiff_p=20.mat');
all_anim_bayes = all_anim_bayes.bf;

all_cat_bayes = load(outdir + 'all_bf_catdiff_p=20.mat');
all_cat_bayes = all_cat_bayes.bf;

aspcat_asp_bayes = load(outdir + 'aspcat_bf_aspdiff_p=20.mat');
aspcat_asp_bayes = aspcat_asp_bayes.bf;

aspcat_cat_bayes = load(outdir + 'aspcat_bf_catdiff_p=20.mat');
aspcat_cat_bayes = aspcat_cat_bayes.bf;

% get timing info
rootdir = "/Users/f004p7b/Documents/australia/objects_project/";
decdir = rootdir + "analysis_bids/data/derivatives/decoding/";
orig_dec = load(decdir + "short_origstim-decoding-allthrough22");
timing = orig_dec.res_animacy.a.fdim.values{1,1};

% put everything together 

fin_all_asp = [];
fin_all_asp(1,:) = mean_all_asp;
fin_all_asp(2,:) = all_asp_bayes;
fin_all_asp(3,:) = timing;
finfile = [outdir + 'diffsbf_all_asp_finalfile.csv'];
writematrix(fin_all_asp,finfile);


fin_all_anim = [];
fin_all_anim(1,:) = mean_all_anim;
fin_all_anim(2,:) = all_anim_bayes;
fin_all_anim(3,:) = timing;
finfile = [outdir + 'diffsbf_all_anim_finalfile.csv'];
writematrix(fin_all_anim,finfile);

fin_all_cat = [];
fin_all_cat(1,:) = mean_all_cat;
fin_all_cat(2,:) = all_cat_bayes;
fin_all_cat(3,:) = timing;
finfile = [outdir + 'diffsbf_all_cat_finalfile.csv'];
writematrix(fin_all_cat,finfile);

fin_aspcat_asp = [];
fin_aspcat_asp(1,:) = mean_aspcat_asp;
fin_aspcat_asp(2,:) = aspcat_asp_bayes;
fin_aspcat_asp(3,:) = timing;
finfile = [outdir + 'diffsbf_aspcat_asp_finalfile.csv'];
writematrix(fin_aspcat_asp,finfile);


fin_aspcat_cat = [];
fin_aspcat_cat(1,:) = mean_aspcat_cat;
fin_aspcat_cat(2,:) = aspcat_cat_bayes;
fin_aspcat_cat(3,:) = timing;
finfile = [outdir + 'diffsbf_aspcat_cat_finalfile.csv'];
writematrix(fin_aspcat_cat,finfile);







