%% temporal generalization analysis - objects australia project
% Nov 1, 2023 AK 


%% set up

close all
clear all
addpath('~/CoSMoMVPA/mvpa/')

tic

%% change this for each participant

for pnum = 1:22

disp([' :) Running temp gen for participant ' num2str(pnum)]) 

%% directories and files 
%this should be the only area that needs to be adjusted when running on other computers 

datdir = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/"; % where data is (bids format)
tempgendir = datdir + "tempgen/"; %where to put temp gen analysis results
T = readtable('/Users/f004p7b/Documents/australia/objects_project/masterlist_YodB.csv'); %reads in master list that has variable information in it
T = sortrows(T,"Var1"); %sorts in alphabetic order (variable 1)

if pnum < 10
    datfile = datdir + "cosmomvpa/sub-0" + pnum + "_cosmomvpa.mat" ;
    svfile = tempgendir + "sub-0" + pnum ;
else
    datfile = datdir + "cosmomvpa/sub-" + pnum + "_cosmomvpa.mat" ;
    svfile = tempgendir + "sub-" + pnum ;
end 

% make the temporal generalization folder location if it doesn't exist
 if ~exist("tempgendir")
        mkdir(tempgendir)
end

%% manipulate data (slicy-slice)

load(datfile)
ds = cosmo_slice(ds,~ds.sa.istarget); %removes the target trials (squares and triangles)
ds.sa.binary = contains(ds.sa.stim,'binarized'); %creates sample attribute that specifies if the trial was using binarized stimuli or not (0 = original, 1 = binarized)
ds.sa.stimname2 = T.Var1(ds.sa.stimnumber); %creates sample attribute that specifies the identity of the trial stimulus
ds.sa.categoryname = T.Var2(ds.sa.stimnumber); %creates sample attribute that specifies the category of each trial stimulus
[~,~,ds.sa.categorynumber] = unique(ds.sa.categoryname); %changes the category names into numbers (1 = bodies, 2 = faces, 3 = manmade object, 4 = natural object)
ds.sa.aspectname = T.Var3(ds.sa.stimnumber);%creates sample attribute that specifies the aspect ratio of each trial stimulus
[~,~,ds.sa.aspectnr] = unique(ds.sa.aspectname); %changes stubby/spiky to number (1 = median, 2 = spikey, 3 = stubby)
ds.sa.conditionnumber = T.Var4(ds.sa.stimnumber); % number of the condition within each category (# between 1-13)
ds.sa.animname = T.Var5(ds.sa.stimnumber); % stimulus animacy status
ds.sa.animacy = double(contains(ds.sa.animname, 'inanimate')); % turns animacy names into numbers (0 = animate, 1 = inanimate)

% need to find the folds for different runs

%% decoding
% goal is to train on each time point and test on all other time points 
% will need it to be aspect ratio once & category once (maybe animacy??)
% also need to do it on binarized and original stimuli separately 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aspect ratio
res_aspr_orig = {}; %where output will go for original stimuli
res_aspr_bin = {}; %where output will go for binarized stimuli
allaspr_out = {}; %try to save both original and binarized stimuli results here
allar_file = svfile + "_asprat_all.mat";


 for b = 0:1
        dslice = cosmo_slice(ds,ds.sa.binary==b & ds.sa.aspectnr>1); %slice the data (based on binarized stimuli and get rid of the medians!)
            if b == 0
                disp('starting aspect ratio temporal decoding for non-binarized stimuli!')
                output_res = res_aspr_orig;
                filename = svfile +  "_asprat_tempgen_orig.mat";
            else
                disp('starting aspect ratio temporal decoding for binarized stimuli!')
                output_res = res_aspr_bin;
                filename = svfile + "_asprat_tempgen_bin.mat";
            end
   
        dslice.sa.targets = dslice.sa.aspectnr ;    
        uc = unique(dslice.sa.sequencenumber); %get the unique names of sequences (runs - should be 20) - will need for folds for training and testing
            for folds = 1:length(uc)
                 disp(['fold #' num2str(folds) ' out of ' num2str(length(uc))])
                 dslice.sa.chunks = ones(length(dslice.sa.sequencenumber),1); %train
                 dslice.sa.chunks(dslice.sa.sequencenumber == uc(folds)) = 2; %test
                 opt.measure=@cosmo_crossvalidation_measure;
                 opt.classifier=@cosmo_classify_lda;
                 opt.dimension='time';
                 opt.output='accuracy';
                 opt.nproc=4;
                 dslice_time=cosmo_dim_transpose(dslice,'time',1);
                 res_apr = cosmo_dim_generalization_measure(dslice_time,opt);
                 output_res{folds} = res_apr;


            end
                save(filename,'output_res')  % don't forget to save next time!! 
                t = b+1;
                allaspr_out{t} = output_res;

 end

                 save(allar_file,'allaspr_out') % will save everything

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% category
res_cat_orig = {}; %where output will go for original stimuli
res_cat_bin = {}; %where output will go for binarized stimuli
allcat_out = {}; %try to save both original and binarized stimuli results here
allcat_file = svfile + "tempgen_cat_all.mat";


 for b = 0:1
        dslice = cosmo_slice(ds,ds.sa.binary==b); %slice the data (based on binarized stimuli)
            if b == 0
                disp('starting category temporal decoding for non-binarized stimuli!')
                output_res = res_cat_orig;
                filename = svfile +  "_cat_tempgen_orig.mat";
            else
                disp('starting category temporal decoding for binarized stimuli!')
                output_res = res_cat_bin;
                filename = svfile + "_cat_tempgen_bin.mat";
            end
   
        dslice.sa.targets = dslice.sa.categorynumber ;    
        uc = unique(dslice.sa.sequencenumber); %get the unique names of sequences (runs - should be 20) - will need for folds for training and testing
            for folds = 1:length(uc)
                 disp(['fold #' num2str(folds) ' out of ' num2str(length(uc))])
                 dslice.sa.chunks = ones(length(dslice.sa.sequencenumber),1); %train
                 dslice.sa.chunks(dslice.sa.sequencenumber == uc(folds)) = 2; %test
                 opt.measure=@cosmo_crossvalidation_measure;
                 opt.classifier=@cosmo_classify_lda;
                 opt.dimension='time';
                 opt.output='accuracy';
                 opt.nproc=4;
                 dslice_time=cosmo_dim_transpose(dslice,'time',1);
                 res_cat = cosmo_dim_generalization_measure(dslice_time,opt);
                 output_res{folds} = res_cat;


            end
                save(filename,'output_res')  
                t = b+1;
                allcat_out{t} = output_res;

 end

                 save(allcat_file,'allcat_out') % will hopefully save everything

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% animacy
res_anim_orig = {}; %where output will go for original stimuli
res_anim_bin = {}; %where output will go for binarized stimuli
allanim_out = {}; %try to save both original and binarized stimuli results here
allanim_file = svfile + "tempgen_anim_all.mat";


 for b = 0:1
        dslice = cosmo_slice(ds,ds.sa.binary==b); %slice the data (based on binarized stimuli)
            if b == 0
                disp('starting animacy temporal decoding for non-binarized stimuli!')
                output_res = res_anim_orig;
                filename = svfile +  "_anim_tempgen_orig.mat";
            else
                disp('starting animacy temporal decoding for binarized stimuli!')
                output_res = res_anim_bin;
                filename = svfile + "_anim_tempgen_bin.mat";
            end
   
        dslice.sa.targets = dslice.sa.animacy ;    
        uc = unique(dslice.sa.sequencenumber); %get the unique names of sequences (runs - should be 20) - will need for folds for training and testing
            for folds = 1:length(uc)
                 disp(['fold #' num2str(folds) ' out of ' num2str(length(uc))])
                 dslice.sa.chunks = ones(length(dslice.sa.sequencenumber),1); %train
                 dslice.sa.chunks(dslice.sa.sequencenumber == uc(folds)) = 2; %test
                 opt.measure=@cosmo_crossvalidation_measure;
                 opt.classifier=@cosmo_classify_lda;
                 opt.dimension='time';
                 opt.output='accuracy';
                 opt.nproc=4;
                 dslice_time=cosmo_dim_transpose(dslice,'time',1);
                 res_anim = cosmo_dim_generalization_measure(dslice_time,opt);
                 output_res{folds} = res_anim;


            end
                save(filename,'output_res')  % don't forget to save next time!! 
                t = b+1;
                allanim_out{t} = output_res;

 end

                 save(allanim_file,'allanim_out') % will hopefully save everything


end 
toc


