%% bayes for difference in temp gen values between binarized and intact stimuli
% AK - 10/17/24

%% Setup
close all;
clearvars;

rootdir = "/Users/f004p7b/Documents/australia/objects_project/";
tempdir = rootdir + "analysis_bids/data/derivatives/tempgen/";
codedir = rootdir + "code/";
addpath(genpath(rootdir + 'BFF_repo'))
addpath(genpath("private/var/folders/0r/3b_8szpj7dq8tjwhbsdmsf9c0000gn/T/Rtmpz35GJw/downloaded_packages/BayesFactor/"))


% load the stacked temporal generalization files 
tmpcat_bin_file = load(tempdir + "all_p=20_bin_cat");
cat_bin_file = tmpcat_bin_file.catbin_all;

tmpcat_orig_file = load(tempdir + "all_p=20_orig_cat");
cat_orig_file = tmpcat_orig_file.catorig_all;

tmpanim_bin_file = load(tempdir + "all_p=20_bin_anim");
anim_bin_file = tmpanim_bin_file.animbin_all;

tmpanim_orig_file = load(tempdir + "all_p=20_orig_anim");
anim_orig_file = tmpanim_orig_file.animorig_all;

tmpasp_bin_file = load(tempdir + "all_p=20_bin_asprat");
asp_bin_file = tmpasp_bin_file.arbin_all;

tmpasp_orig_file = load(tempdir + "all_p=20_orig_asprat");
asp_orig_file = tmpasp_orig_file.arorig_all;

% calculate the differences 
cat_diff_all{20,1} = {};
asp_diff_all{20,1} = {};
anim_diff_all{20,1} = {};

for part = 1:size(anim_bin_file,1)
    % isolate each participant's temp gen results & get rid of extra
    % dimension
    anim_int = squeeze(anim_orig_file{part});
    anim_bin = squeeze(anim_bin_file{part}); 
    asp_int = squeeze(asp_orig_file{part});
    asp_bin = squeeze(asp_bin_file{part}); 
    cat_int = squeeze(cat_orig_file{part});
    cat_bin = squeeze(cat_bin_file{part}); 

    %set temp. mat for each participant 
    anim_diff_temp = zeros(181, 181);
    asp_diff_temp = zeros(181, 181);
    cat_diff_temp = zeros(181, 181); 

    %have to calculate the difference at each time point in each matrix
    for tim1 = 1:size(cat_bin,2)
        for tim2 = 1:size(cat_bin,2)
            val_cat_int = cat_int(tim1,tim2);   %category intact decoding value
            val_cat_bin = cat_bin(tim1,tim2);   %category binarized decoding value
            val_asp_int = asp_int(tim1,tim2);  %aspect ratio intact decoding value
            val_asp_bin = asp_bin(tim1,tim2);   %aspect ratio binarized decoding value
            val_anim_int = anim_int(tim1,tim2); %animacy intact decoding value
            val_anim_bin = anim_bin(tim1,tim2); %animacy binarized decoding value

            val_cat_diff_tmp = val_cat_int - val_cat_bin;
            val_asp_diff_tmp = val_asp_int - val_asp_bin;
            val_anim_diff_tmp = val_anim_int - val_anim_bin; 

            anim_diff_temp(tim1,tim2) = val_anim_diff_tmp;
            asp_diff_temp(tim1,tim2) = val_asp_diff_tmp;
            cat_diff_temp(tim1,tim2) = val_cat_diff_tmp; 
        end
    end

 cat_diff_all{part, 1} = cat_diff_temp; 
 anim_diff_all{part, 1} = anim_diff_temp; 
 asp_diff_all{part, 1} = asp_diff_temp;

end

% save for bayes
anim_outfile = [tempdir + 'animtempgen_diffs_p=20.mat']; 
asp_outfile = [tempdir + 'asprattempgen_diffs_p=20.mat']; 
cat_outfile = [tempdir + 'cattempgen_diffs_p=20.mat']; 

save(anim_outfile,"anim_diff_all")
save(asp_outfile,"asp_diff_all")
save(cat_outfile,"cat_diff_all")

%% bayes analysis 
% need a file that is 20 x 181 x 181
mval = ['mu=0, rscale="medium",nullInterval=c(-0.5,0.5)']; % set parameters 

for dims = 1:3
    if dims == 1
        % category
            %file = cell2mat(cat_diff_all); % make this into a matrix so that we can do math - converts 20 x 1 cell (each cell is 181 x 181) to a 3620 x 181 matrix
        test1 = cell2table(cat_diff_all);
        test2 = table2struct(test1);
        file = [];
        for ugh = 1:size(test2,1)
            file(ugh,:,:) = test2(ugh).cat_diff_all;
        end
        outfile = [tempdir + 'bf_tempgen_catdiff_p=20.mat']; 
        disp('starting bayes for cat differences')
    elseif dims == 2
        % aspect ratio 
        test1 = cell2table(asp_diff_all); 
        test2 = table2struct(test1);
        file = [];
        for ugh = 1:size(test2,1)
            file(ugh,:,:) = test2(ugh).asp_diff_all;
        end
        outfile = [tempdir + 'bf_tempgen_aspratdiff_p=20.mat']; 
        disp('starting bayes for aspect ratio differences')
    elseif dims == 3
        % animacy
        test1 = cell2table(anim_diff_all);
        test2 = table2struct(test1);
        file = [];
        for ugh = 1:size(test2,1)
            file(ugh,:,:) = test2(ugh).anim_diff_all;
        end
        outfile = [tempdir + 'bf_tempgen_animdiff_p=20.mat']; 
        disp('starting bayes for animacy differences')
    end

%mat_size = size(file,2);
 bf = zeros(size(file,[2,3])); %make a blank matrix that is the size of 181 x 181 where bayes info will be stored

    for bf_diff = 1:size(bf,2)
        disp(['running bayes for timepoint ' num2str(bf_diff), ' out of ' num2str(size(bf,2))])
        tic
        X = file(:,:,bf_diff)'; %changed this because first time this ran, it saved a 1 x 181 instead of 181 x 181
        bf(bf_diff,:) = bayesfactor_R_wrapper(X,'args',mval,'returnindex',2);
    end

    save(outfile,"bf")
end
toc
    







