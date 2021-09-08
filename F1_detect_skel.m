%Code used to detect filaments in image files
%Change the variables under common variables to optimize detection 
%depending on the image
%F1 and F2 contain codes adapted from @cwchung1996/AFMFibrilPolyMorph

%common variables to change:
new = 1;      %set to 1 to load new tiff image file, 0 to keep current image
thres =0.4;   %set threshold using fraction of max value (default of 0.25 = 25%)
gus = 4;      %level of gaussian filtering for detection only

%
top_hat = 0;  %leveling option 
hd=0;         %high density? set to 1 to reduce branch numbers
mbl = 50;     %min branch length used to reduce branch length

%%
if new ==1
clearvars -except mbl hd top_hat thres new gus
[f,path] = uigetfile('*.tif');
f = fullfile(path,f); %tif (in nm) filename
info = imfinfo(f);  n = numel(info);
A = imread(f);
else
end
load('lutafm.mat')
%%
sz = size(A);
ndata =1;
count=1;
if (top_hat~=0)
 se = strel('disk',top_hat);
 Ap = imtophat(A,se);
else
    Ap = A; 
end

if (gus~=0)
    Ag = imgaussfilt(Ap,gus);
    else 
    Ag = Ap;    
    end
 At=(Ag>(thres*max(Ag(:))));

       [clusters, nclus]=bwlabeln(At);
        c1=1;c2=length(clusters);
         skelimg=zeros(size(A));
         

        s1=bwmorph(clusters, 'thin', Inf) ;
        s1 = bwskel(logical(s1),'MinBranchLength',mbl);
        
s1 = bwmorph(s1, 'spur',2) ;
if (hd==1)
s2 = bwmorph(s1, 'branchpoints') ;
pos = s2>0;
s1(pos) = 0;
else if (hd==2)
        s2 = bwmorph(s1, 'branchpoints') ;
        s2 = bwmorph(s2, 'dilate',1) ;
        pos = s2>0;
        s1(pos) = 0;  
    else
    end
end
    s1 = bwmorph(s1,'clean');
[clusters, nclus]=bwlabeln(s1);
figure('Position',[10 1000 1500 400])
tiledlayout(1,3, 'Padding', 'none', 'TileSpacing', 'none');  
nexttile  
imagesc(Ap)
nexttile 
imagesc(At)
%colormap(lutafm)
nexttile 
imagesc(clusters)
