%% decoding script for the australia objects project 
% using the output from the preprocessing scripts with shortened epochs (.8
% seconds after stimulus onset instead of 1.2 seconds after stimulus onset)
    % Nov 8 2023 - AK

clear all

%% participant number - change each time you run

pnum = 22;


%% files and set up 

if pnum < 10 
    matsource = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/cosmomvpa/sub-0" +pnum + "_cosmomvpa.mat" ; % output of the preprocessing script
else
    matsource = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/cosmomvpa/sub-" +pnum + "_cosmomvpa.mat" ; % output of the preprocessing script
end

load(matsource) 
ds = cosmo_slice(ds,~ds.sa.istarget); % takes out "target" trials when responding to triangle or rectangle -  resulting ds samples = 8320 (stimuli events/presented) x 16,704 (channels {64} x timepoints {261})
T = readtable('/Users/f004p7b/Documents/australia/objects_project/masterlist_YodB.csv'); %reads in master list that has variable information in it
T = sortrows(T,"Var1"); %sorts in alphabetic order (variable 1)
ds.sa.binary = contains(ds.sa.stim,'binarized'); %creates sample attribute that specifies if the trial was using binarized stimuli or not (0 = original, 1 = binarized)
ds.sa.stimname2 = T.Var1(ds.sa.stimnumber); %creates sample attribute that specifies the identity of the trial stimulus
ds.sa.categoryname = T.Var2(ds.sa.stimnumber); %creates sample attribute that specifies the category of each trial stimulus
[~,~,ds.sa.categorynumber] = unique(ds.sa.categoryname); %changes the category names into numbers (1 = bodies, 2 = faces, 3 = manmade object, 4 = natural object)
ds.sa.aspectname = T.Var3(ds.sa.stimnumber);%creates sample attribute that specifies the aspect ratio of each trial stimulus
[~,~,ds.sa.aspectnr] = unique(ds.sa.aspectname); %changes stubby/spiky to number (1 = median, 2 = spikey, 3 = stubby)
ds.sa.conditionnumber = T.Var4(ds.sa.stimnumber); % number of the condition within each category (# between 1-13)
ds.sa.animname = T.Var5(ds.sa.stimnumber); % stimulus animacy status
ds.sa.animacy = double(contains(ds.sa.animname, 'inanimate')); % turns animacy names into numbers (0 = animate, 1 = inanimate)


%% image (stimulus identity decoding)
ds.sa.targets = ds.sa.stimnumber; %specifies targets as the stimulus number (what you are decoding)
ds.sa.chunks = ds.sa.yokedsequencenumber; %chunks that are used are the yoked sequence number (the sequence number that is common across all participants)
nh = cosmo_interval_neighborhood(ds,'time','radius',0); %specifies the neighborhood for decoding searchlight (dataset to use, dimension to compute the neighborhood (time), radius = 0 means not binning time)      


% as a note, in res_cell_image, there are 2 cells - the first cell corresponds to nonbinarized images and the second to binarized images

res_cell_image ={}; %where the decoding results will be 
for b=0:1 %for both binarized = 0 and 1
    dsb = cosmo_slice(ds,ds.sa.binary==b); %what is used to break up the dataset - in this case, if stimulus is binarized vs original
    dsb = cosmo_average_samples(dsb,'count',4); % averages subsets of the dataset by 4 samples per average
    ma = {};
    ma.classifier = @cosmo_classify_lda; %type of classifier being used (linear discriminant analysis classifier)
    ma.partitions = cosmo_nfold_partitioner(dsb); %creates the folds/partition scheme for the training/testing of the decoder - uses dsb.sa.chunks (which is 1 x 1040 - 1040 x 8 = total number of trials {potentially the 4 specified above x 2 because there is binarized and nonbinarized})
    ma.nproc=4;
    res = cosmo_searchlight(dsb,nh,@cosmo_crossvalidation_measure,ma); %actually performing the decoding
    res_cell_image{b+1} = res; % which cell are the results put in original = first cell :) 
end


%% category - see general comments above for the commands below
ds.sa.targets = ds.sa.categorynumber;  %specifies targets as the category number (what you are decoding)
ds.sa.chunks = ds.sa.yokedsequencenumber;
nh = cosmo_interval_neighborhood(ds,'time','radius',0);
res_cell_cat ={};
for b=0:1
    dsb = cosmo_slice(ds,ds.sa.binary==b);
    ma = {};
    ma.classifier = @cosmo_classify_lda;
    ma.partitions = cosmo_nfold_partitioner(dsb);
    ma.nproc=4;
    res = cosmo_searchlight(dsb,nh,@cosmo_crossvalidation_measure,ma);
    res_cell_cat{b+1} = res;
end

%% aspect
ds.sa.targets = ds.sa.aspectnr;
ds.sa.chunks = ds.sa.yokedsequencenumber;
nh = cosmo_interval_neighborhood(ds,'time','radius',0);
res_cell_asp ={};
for b=0:1
    dsb = cosmo_slice(ds,ds.sa.binary==b & ds.sa.aspectnr>1); %remove medians
    %dsb = cosmo_average_samples(dsb,'count',4);
    ma = {};
    ma.classifier = @cosmo_classify_lda;
    ma.partitions = cosmo_nfold_partitioner(dsb);
    ma.nproc=4;
    res = cosmo_searchlight(dsb,nh,@cosmo_crossvalidation_measure,ma);
    res_cell_asp{b+1} = res;
end

%% animacy

ds.sa.targets = ds.sa.animacy; %if the stimulus is animate vs inanimate
ds.sa.chunks = ds.sa.yokedsequencenumber;
nh = cosmo_interval_neighborhood(ds,'time','radius',0);
res_cell_anim ={};
for b=0:1
    dsb = cosmo_slice(ds,ds.sa.binary==b); 
    ma = {};
    ma.classifier = @cosmo_classify_lda;
    ma.partitions = cosmo_nfold_partitioner(dsb);
    ma.nproc=4;
    res = cosmo_searchlight(dsb,nh,@cosmo_crossvalidation_measure,ma);
    res_cell_anim{b+1} = res;
end

%% save file

if pnum < 10
    decoding_files = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/decoding/sub-0" + pnum + "_res_decoding.mat";
    save(decoding_files,'res_cell_asp','res_cell_cat','res_cell_image', 'res_cell_anim')
else
    decoding_files = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/decoding/sub-" + pnum + "_res_decoding.mat";
    save(decoding_files,'res_cell_asp','res_cell_cat','res_cell_image', 'res_cell_anim')
end

%% plot
figure(1);clf
subplot(4,1,1)
tv = res_cell_image{1}.a.fdim.values{1};
for i=1:2
    plot(tv,res_cell_image{i}.samples);hold on
end
plot(tv,1/52+0*tv,'k--')
legend({'original','binarized','chance'})
xlim(minmax(tv))
title('decoding image (52 images)')

subplot(4,1,2)
tv = res_cell_cat{1}.a.fdim.values{1};
for i=1:2
    plot(tv,res_cell_cat{i}.samples);hold on
end
plot(tv,1/4+0*tv,'k--')
legend({'original','binarized','chance'})
xlim(minmax(tv))
title('decoding category (face/body/manmade/natural)')

subplot(4,1,3)
tv = res_cell_asp{1}.a.fdim.values{1};
for i=1:2
    plot(tv,res_cell_asp{i}.samples);hold on
end
plot(tv,1/2+0*tv,'k--')
legend({'original','binarized','chance'})
xlim(minmax(tv))
title('decoding aspect (stubby/spikey)')

subplot(4,1,4)
tv = res_cell_anim{1}.a.fdim.values{1};
for i=1:2
    plot(tv,res_cell_anim{i}.samples);hold on
end
plot(tv,1/2+0*tv,'k--')
legend({'original','binarized','chance'})
xlim(minmax(tv))
title('decoding animacy (animate/inanimate)')

%% save plot
fn = '/Users/f004p7b/Documents/australia/objects_project/plots/decoding_figures/';
if pnum < 10
    tempname = "sub-0" + pnum + "_decoding.png";
else
    tempname = "sub-" + pnum + "_decoding.png";
end

tn = fn + tempname;

print(gcf,'-dpng','-r500',tn)
im=imread([tn]);
[i,j]=find(mean(im,3)<255);margin=0;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[tn '.png'],'png');
 % note: imwrite sometimes throws a weird error, but it still saves and
 % everything is fine so ignore it 




