%% bayes for decoding differences between intact and binarized stimuli 
% AK 10/17/24

% also added calculating standard deviation on 10/24/24

%% Setup
close all;
clearvars;

rootdir = "/Users/f004p7b/Documents/australia/objects_project/";
decdir = rootdir + "analysis_bids/data/derivatives/decoding/";
codedir = rootdir + "code/";
addpath(genpath(rootdir + 'BFF_repo'))
addpath(genpath("private/var/folders/0r/3b_8szpj7dq8tjwhbsdmsf9c0000gn/T/Rtmpz35GJw/downloaded_packages/BayesFactor/"))

% load stacked decoding results (decoding across all participants, not averaged)
orig_dec = load(decdir + "short_origstim-decoding-allthrough22");
bin_dec = load(decdir + "short_binstim-decoding-allthrough22");

% separate results from big files 
%image_orig_file = orig_dec.res_image.samples;
%image_bin_file = bin_dec.res_bin_image.samples;
cat_orig_file = orig_dec.res_cat.samples;
cat_bin_file = bin_dec.res_bin_cat.samples;
anim_orig_file = orig_dec.res_animacy.samples;
anim_bin_file = bin_dec.res_bin_animacy.samples;
asp_orig_file = orig_dec.res_asp.samples;
asp_bin_file = bin_dec.res_bin_asp.samples;


%% standard deviation and standard error 

%sd and save
cat_orig_sd = std(cat_orig_file);
cat_bin_sd = std(cat_bin_file);
anim_orig_sd = std(anim_orig_file);
anim_bin_sd = std(anim_bin_file);
asp_orig_sd = std(asp_orig_file);
asp_bin_sd = std(asp_bin_file);

sd_cat_origfile = [decdir + 'catdecoding_sd_orig.csv'];
sd_cat_binfile = [decdir + 'catdecoding_sd_bin.csv'];
sd_anim_origfile = [decdir + 'animdecoding_sd_orig.csv'];
sd_anim_binfile = [decdir + 'animdecoding_sd_bin.csv'];
sd_asp_origfile = [decdir + 'aspdecoding_sd_orig.csv'];
sd_asp_binfile = [decdir + 'aspdecoding_sd_bin.csv'];

writematrix(cat_orig_sd,sd_cat_origfile);
writematrix(cat_bin_sd,sd_cat_binfile);
writematrix(asp_orig_sd,sd_asp_origfile);
writematrix(asp_bin_sd,sd_asp_binfile);
writematrix(anim_orig_sd,sd_anim_origfile);
writematrix(anim_bin_sd,sd_anim_binfile);

%se and save
cat_orig_se = cat_orig_sd/sqrt(size(cat_orig_file,1));
cat_bin_se = cat_bin_sd/sqrt(size(cat_orig_file,1));
anim_orig_se = anim_orig_sd/sqrt(size(cat_orig_file,1));
anim_bin_se = anim_bin_sd/sqrt(size(cat_orig_file,1));
asp_orig_se = asp_orig_sd/sqrt(size(cat_orig_file,1));
asp_bin_se = asp_bin_sd/sqrt(size(cat_orig_file,1));

se_cat_origfile = [decdir + 'catdecoding_se_orig.csv'];
se_cat_binfile = [decdir + 'catdecoding_se_bin.csv'];
se_anim_origfile = [decdir + 'animdecoding_se_orig.csv'];
se_anim_binfile = [decdir + 'animdecoding_se_bin.csv'];
se_asp_origfile = [decdir + 'aspdecoding_se_orig.csv'];
se_asp_binfile = [decdir + 'aspdecoding_se_bin.csv'];

writematrix(cat_orig_se,se_cat_origfile);
writematrix(cat_bin_se,se_cat_binfile);
writematrix(anim_orig_se,se_anim_origfile );
writematrix(anim_bin_se,se_anim_binfile );
writematrix(asp_orig_se,se_asp_origfile);
writematrix(asp_bin_se,se_asp_binfile);

%% get difference score for each participant at each time point for each dimension

cat_diff_all = [];
asp_diff_all = [];
anim_diff_all = [];

for part = 1:size(cat_bin_file,1)
    for tim = 1:size(cat_bin_file,2)
        cat_int = cat_orig_file(part,tim);  %category intact decoding value
        cat_bin = cat_bin_file(part,tim);   %category binarized decoding value
        asp_int = asp_orig_file(part,tim);  %aspect ratio intact decoding value
        asp_bin = asp_bin_file(part,tim);   %aspect ratio binarized decoding value
        anim_int = anim_orig_file(part,tim); %animacy intact decoding value
        anim_bin = anim_bin_file(part,tim); %animacy binarized decoding value

        cat_diff_tmp = cat_int - cat_bin;
        asp_diff_tmp = asp_int - asp_bin;
        anim_diff_tmp = anim_int - anim_bin; 

        cat_diff_all(part,tim) = cat_diff_tmp;
        asp_diff_all(part,tim) = asp_diff_tmp;
        anim_diff_all(part,tim) = anim_diff_tmp;

    end
end

% save difference files for plotting purposes

anim_outfile = [decdir + 'new_anim_diffs_p=20.mat']; 
asp_outfile = [decdir + 'new_asprat_diffs_p=20.mat']; 
cat_outfile = [decdir + 'new_cat_diffs_p=20.mat']; 

save(anim_outfile,"anim_diff_all")
save(asp_outfile,"asp_diff_all")
save(cat_outfile,"cat_diff_all")

%% stats
% note:
% (from: https://www.rdocumentation.org/packages/BayesFactor/versions/0.9.12-4.2/topics/ttestBF)
%   in the bayesfactor_R_wrapper function:
%       mu = the null value of the mean, or the mean difference (for
%            one-sample, and paired tests)
%       rscale = prior scale; presets can be given as strings (like
%            "medium")
%       nullInterval = optional vector of length 2 containing lower 
%            and upper bounds of an interval hypothesis to test, in standardized units
%            (this is the interval that Bayes factors are being calculated
%            in, not an interval of time! for example, if nullInterval=c(0.5,Inf) this
%            is telling matlab to calculate bayes factors between the
%            values of 0.5 (which is evidence for the null) to infinity
%            (values that are 1 and higher are evidence for the alternative
%            hypothesis)

mval = ['mu=0, rscale="medium",nullInterval=c(-0.5,0.5)'];

for a = 1:3
    if a==1
        file = cat_diff_all;
        outfile = [decdir + 'new_bf_catdiff_p=20.mat']; 
        disp('starting bayes for cat differences')
    elseif a==2
        file = anim_diff_all;
        outfile = [decdir + 'new_bf_animdiff_p=20.mat']; 
        disp('starting bayes for anim differences')  
    elseif a == 3 
        file = asp_diff_all;
        outfile = [decdir + 'new_bf_aspratdiff_p=20.mat']; 
        disp('starting bayes for asp ratio differences')
    end

    bf = zeros(size(file,[2,3]))'; %make a blank matrix that is the size of 181 x 181 where bayes info will be stored
    
    %for bf_diff = 1:size(bf,2)
        %disp(['running bayes for timepoint ' num2str(bf_diff), ' out of ' num2str(size(bf,2))])
        %tic
        X = file';
        bf(1,:) = bayesfactor_R_wrapper(X,'args',mval,'returnindex',2);
    %end
    
    save(outfile,"bf")

end

%% for python figures - make one file that has decoding differences, bayes,
% time for each dimension

% average differences
avg_cat_diff_all = mean(cat_diff_all);
avg_asp_diff_all = mean(asp_diff_all) ;
avg_anim_diff_all = mean(anim_diff_all);

% bayes
bayes_asprat = load('/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/decoding/new_bf_aspratdiff_p=20.mat');
bf_asp = bayes_asprat.bf;

bayes_anim = load('/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/decoding/new_bf_animdiff_p=20.mat');
bf_anim = bayes_anim.bf;

bayes_cat = load('/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/decoding/new_bf_catdiff_p=20.mat');
bf_cat = bayes_cat.bf;

% time info
timing = orig_dec.res_animacy.a.fdim.values{1,1};

% put it together 
cat_all = [];
cat_finfile = [decdir + 'new_cat_diffsbf_finalfile.csv'];
anim_all = [];
anim_finfile = [decdir + 'new_anim_diffsbf_finalfile.csv'];
asp_all = [];
asp_finfile = [decdir + 'new_asp_diffsbf_finalfile.csv'];


cat_all(1:181,1) = avg_cat_diff_all';
cat_all(1:181,2) = bf_cat(1,1:181)';
cat_all(1:181,3) = timing';
% save(cat_finfile,"cat_all");
writematrix(cat_all,cat_finfile)

anim_all(1:181,1) = avg_anim_diff_all';
anim_all(1:181,2) = bf_anim(1,1:181)';
anim_all(1:181,3) = timing';
%save(anim_finfile,"anim_all");
writematrix(anim_all,anim_finfile)

asp_all(1:181,1) = avg_asp_diff_all';
asp_all(1:181,2) = bf_asp(1,1:181)';
asp_all(1:181,3) = timing';
%save(asp_finfile,"asp_all");
writematrix(asp_all,asp_finfile)


