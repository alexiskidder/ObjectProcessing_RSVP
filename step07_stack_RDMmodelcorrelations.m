%% stacking model RDM results and averaging together
% November 15, 2023 - AK


%% set up and load 

clear all

basedir = "/Users/f004p7b/Documents/australia/objects_project/";
datadir = basedir + "analysis_bids/data/derivatives/modelrdm/" ; %where RDM model correlation outputs are located
figdir = basedir + "modelrdmcorr_figures/"; % where you want the figure to go

files = dir(datadir + "*.mat"); %where all the preprocessed data is (get names of the files)
svfile = datadir + "averageRDMmodelcorr_p=" + length(files) + ".mat"; % name of file you will save data to
plotfile = figdir + "averageRDMmodelcorr_p=" + length(files) + ".png"; %name of plot file

all_origcorr = []; % stack non-binarized correlations 
all_bincorr = [] ; % stack binarized correlations

%% get all participant files 

for k = 1:length(files)

    basefil = files(k).name ; %this is the name of the file you will need 
    pnum = basefil(5:6); %gets the participant number you are working on from the file you are loading
    corr_file = datadir + basefil;
    pfile = load(corr_file);
    pfile = pfile.all_corr;
    all_origcorr(k,:,:) = pfile{1};
    all_bincorr(k,:,:) = pfile{2};
end

%% average

reframe_orig = permute(all_origcorr,[2 3 1]); % want participant # to be 3rd dimension
res_origcorr = mean(reframe_orig, 3); % average across 3rd dimension (participant)
reframe_bin = permute(all_bincorr,[2 3 1]);
res_bincorr = mean(reframe_bin,3);

%% save

all_average{1} = res_origcorr;
all_average{2} = res_bincorr;

save(svfile,"all_average",'-mat');

%% plot 

figure(1); clf
 % original
   subplot(2,1,1)
   plot(all_average{1}); hold on
   yline(0,'--k')
   xlim([0 181])
   ylim([-.2 .4])
   legend({'aspect ratio','animacy','category'})
   title('model fits - nonbinarized stimuli')

 % binarized
   subplot(2,1,2)
   plot(all_average{2}); hold on
   yline(0,'--k')
   xlim([0 181])
   ylim([-.2 .4])
   legend({'aspect ratio','animacy','category'})
   title('model fits - binarized stimuli')

 % save

Image = getframe(gcf) ; % gets info from the opened figure
imwrite(Image.cdata, plotfile)














