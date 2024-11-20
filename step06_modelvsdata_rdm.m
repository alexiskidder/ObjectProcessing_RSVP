%% updated modelvsdata_rdm  
% nov 9 2023 - ak

clear all 

%% set up and load data

% directories 
datadir = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/rdm/" ;% where data is from
modeldir = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/modelRDMS/" ; % where model RDMs are located
outdir = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/modelrdm/" ; %where output will be placed

% list of participants
parts = dir(datadir + "*.mat"); %where all the preprocessed data is (get names of the files)

% model files 
models = dir(modeldir + "*.mat"); %get the available model names

animacyfile = modeldir + models(1).name ;
aspratfile = modeldir + models(2).name ; 
catfile = modeldir + models(3).name ; 

% manipulate models
    % Note: can't do a model for image identity because it would just all be
    % different (unless potentially get behavioral data that indicates which
    % objects physically look more similar to each other and therefore are more
    % likely to be confused easier?)


%% loop-de-loop

for k = 20%1:length(parts)

    basefil = parts(k).name ; %this is the name of the file you will need  
    pnum = basefil(5:6); %gets the participant number you are working on from the file you are loading
    disp("starting participant " + pnum)

    % set and load files
    svfile = outdir + "sub-" + pnum + "_modelRDM_correlations.mat";
    plotfile = "sub-" + pnum + "_modelRDM_correlationplot.png";
    rdms_all = load(datadir + basefil);
    origrdm = rdms_all.rdm_store{1}; %non-binarized stimuli
    binrdm = rdms_all.rdm_store{2}; %binarized stimuli


%% manipulate the data into a more useable format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% empirical data
% Notes: 
    % 1: you don't want to use data along the diagnoal! 

  % nonbinarized stimuli 
    uporig = permute(origrdm,[3 1 2]); %resulting size should be [181 52 52] - manipulate the matrices 
    loweridx = find(tril(ones(size(origrdm(:,:,1))),-1)); %gets the position of comparisons that you want in the matrix above 
    uporig = uporig(:,loweridx); %should result in 181 x 1326 (timepoints x data comparisons)

 % binarized stimuli
    bin_uporig = permute(binrdm,[3 1 2]); %resulting size should be [181 52 52] - manipulate the matrices 
    bin_loweridx = find(tril(ones(size(binrdm(:,:,1))),-1)); %gets the position of comparisons that you want in the matrix above 
    bin_uporig = bin_uporig(:,bin_loweridx); %should result in 181 x 1326 (timepoints x data comparisons)
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% models

allmodels = [] ; % where all models will be placed (should end up as 1326 x 3)

 % aspect ratio
asmodel = load(aspratfile);
armodel = asmodel.binar_model;
armodel = armodel(loweridx); 
allmodels(1:1326,1) = armodel;

 % animacy
animmodel = load(animacyfile);
anmodel = animmodel.binanimat_model;
animodel = anmodel(loweridx);
allmodels(1:1326,2) = animodel;

 % category 
catmodel = load(catfile);
cmodel = catmodel.bincat_model;
camodel = cmodel(loweridx);
allmodels(1:1326,3) = camodel;


%% actually do the analysis now! 

% Note:     
    % 1: use a spearman (rank correlation, because not assuming there is a linear
    % relationship)
    % 2: output = correlation at each time point for each model (row 1 =
    % aspect ratio, row 2 = animacy, row 3 = category)

% binarized

disp("calcualting correlation between models and binarized stimuli")
bin_corr = corr(bin_uporig',allmodels,'type','Spearman');
bin_corr_t = bin_corr'; % will result in 3 x 181 (model correlation x timepoint)

% original 

disp("calcualting correlation between models and original stimuli")
orig_corr = corr(uporig',allmodels,'type','Spearman');
orig_corr_t = orig_corr'; % will result in 3 x 181 (model correlation x timepoint)

% save
all_corr{1} = orig_corr;
all_corr{2} = bin_corr;

save(svfile,"all_corr",'-mat')

%% plot 

figure(1); clf
 % original
   subplot(2,1,1)
   plot(orig_corr); hold on
   yline(0,'--k')
   xlim([0 181])
   ylim([-.2 .4])
   legend({'aspect ratio','animacy','category'})
   title('model fits - nonbinarized stimuli')

 % binarized
   subplot(2,1,2)
   plot(bin_corr); hold on
   yline(0,'--k')
   xlim([0 181])
   ylim([-.2 .4])
   legend({'aspect ratio','animacy','category'})
   title('model fits - binarized stimuli')


%% save plot
fn = '/Users/f004p7b/Documents/australia/objects_project/modelrdmcorr_figures/';
tn = fn + plotfile;
Image = getframe(gcf) ; % gets info from the opened figure
imwrite(Image.cdata, tn)


 close all

end












