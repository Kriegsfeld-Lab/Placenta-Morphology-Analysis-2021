warning off


Image_Selection = 'How many images do you want to run?';
Image_Range = input(Image_Selection);

for Image_Name = 1:Image_Range
Prompt = 'What is the name of the file exactly? ';      %type the file name excluding the ".tif" tag
File_Name(Image_Name,:) = input(Prompt,'s');    
end

%%
for Script_Cycle = 1:Image_Range
fprintf('Code Started at %s\n', datestr(now,'HH:MM:SS'))
%Creating name for Script Run
Title = File_Name(Script_Cycle,:);
%%
PlacentaDirectory = 'C:\Users\samfi\Desktop\downloads\';          %choose directory where input photo exists
PlacentaName = [PlacentaDirectory Title]; %here
ImageType = '.tif';
PlacentaFile = [PlacentaName ImageType];
he = imread(PlacentaFile);


%%
%Start of code
resized_he = he;
he_conv = rgb2gray(he);
he_conv = imbinarize(he_conv);
he_conv = ~he_conv;
he_conv = bwareaopen(he_conv,10000);
he_conv(:,:,2) = he_conv;
he_conv(:,:,3) = he_conv(:,:,1);
resized_he(he_conv == 0) = 0;



he_resize = imresize(resized_he,0.25);  %Size down the pixel count, remember to size up your results to account for this
Outline = he_resize;
Binary_Outline = rgb2gray(Outline);
Binary_Outline = imbinarize(Binary_Outline);
Binary_Outline = imfill(Binary_Outline,'holes');
Binary_Outline = imerode(Binary_Outline,strel('disk',10));
Binary_Inverse = ~Binary_Outline;
[indConv,mapConv] = rgb2ind(he_resize,16);
he_resize = ind2rgb(indConv,mapConv);
%%
cform = makecform('srgb2cmyk');
lab_he = applycform(he_resize,cform);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);
nColors = 8;

[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','cosine', ...
	'Replicates',1);

pixel_labels = reshape(cluster_idx,nrows,ncols);
sub_clusters = cell(1,8);
rgb_label = repmat(pixel_labels,[1 1 3]);

for k = 1:nColors
	color = he_resize;
	color(rgb_label ~= k) = 0;
	sub_clusters{k} = color;
end

figure
montage({sub_clusters{1},sub_clusters{2},sub_clusters{3},sub_clusters{4},sub_clusters{5},sub_clusters{6},sub_clusters{7},sub_clusters{8}})

%%
MainCluster = 'Which image contains all three colors in regions?';
Cluster_chosen = input(MainCluster);
he_cluster = sub_clusters{Cluster_chosen};

close
%%
cform = makecform('srgb2cmyk');
lab_he = applycform(he_cluster,cform);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);
nColors = 8;

[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','cosine', ...
	'Replicates',1);

pixel_labels = reshape(cluster_idx,nrows,ncols);
main_clusters = cell(1,8);
rgb_label = repmat(pixel_labels,[1 1 3]);

for k = 1:nColors
	color = he_cluster;
	color(rgb_label ~= k) = 0;
	main_clusters{k} = color;
end
figure
montage({main_clusters{1},main_clusters{2},main_clusters{3},main_clusters{4},main_clusters{5},main_clusters{6},main_clusters{7},main_clusters{8}})



%%
LZ_Question = 'Which One is The LZ?';
JZ_Question = 'Which One is The JZ?';
Decidua_Question = 'Which One is The Decidua?';

LZ = input(LZ_Question);
JZ = input(JZ_Question);
Decidua = input(Decidua_Question);

close
%%
Section{1} = main_clusters{LZ};
Section{2} = main_clusters{JZ};
Section{3} = main_clusters{Decidua};

%%
%Grab large empty sections from JZ Zone to add to the mising portions of color
Background = Outline;

for mm = 1:size(Binary_Inverse, 1)
for nn = 1:size(Binary_Inverse, 2)
    if Binary_Inverse(mm,nn,1) > 0
        gsc = (255 + (0.0*Background(mm,nn,1))) + (255 + (0.0*Background(mm,nn,2))) + (255 + (0.0*Background(mm,nn,3)));
        Background(mm,nn,:) = [gsc gsc gsc];
    end
end
end





Background = rgb2gray(Background);
Background = imbinarize(Background, 0.1);
Background = ~Background;

Background_Tissue = cell(1,3);
for n = 1:3
Background_Tissue{n} = Background;
end
%%
%Crop section for LZ out to have more accurate pre processing 
Mask = cell(1,3);
EmptyTissue = cell(1,3);
for n = 1:3
figure,imshow(Section{n})
h = drawfreehand;
Mask{n} = h.createMask;
close
EmptyTissue{n} = Mask{n};
end




for n = 1:3
    EmptyTissue{n} = ~EmptyTissue{n};
    Background_Tissue{n} = Background_Tissue{n} - EmptyTissue{n};
end
%%
%Cut mask out of sections
for n = 1:3
Mask{n}(:,:,2) = Mask{n};
Mask{n}(:,:,3) = Mask{n}(:,:,1);
Section{n}(Mask{n} == 0) = 0;
end

%%
Shift = [25 25; 25 0; 25 -25; 0 -25; -25 -25; -25 0; -25 25; 0 25];

    eight_bit_shift = cell(3,8);
for m = 1:3
    for n = 1:8
    eight_bit_shift{m,n} = circshift(Section{m},Shift(n));
    end
end

Shifted = cell(3,1);
Shifted{1} = Section{1};
Shifted{2} = Section{2};
Shifted{3} = Section{3};

for m = 1:3
    for n = 1:8
    Shifted{m} = imfuse(Shifted{m},eight_bit_shift{m,n},'blend');
    end
end
%%

for n = 1:3

[L, N] = superpixels(Shifted{n},5000);
BW = boundarymask(L);


outputImage = zeros(size(Shifted{n}),'like',Shifted{n});
idx = label2idx(L);
numRows = size(Shifted{n},1);
numCols = size(Shifted{n},2);
for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    outputImage(redIdx) = mean(Shifted{n}(redIdx));
    outputImage(greenIdx) = mean(Shifted{n}(greenIdx));
    outputImage(blueIdx) = mean(Shifted{n}(blueIdx));
end    

Segments{n} = outputImage;
end
%%
%Test
Gray = cell(1,3);
for n = 1:3
Gray{n} = rgb2gray(Section{n});
end

%%
Binary = cell(1,3);
for n = 1:3
Binary{n} = rgb2gray(Segments{n});
Binary{n} = imbinarize(Binary{n},'adaptive','ForegroundPolarity','dark','Sensitivity',0.2);
end

%%
%Fill sections then continue to grow areas with imdilate
Dilated = cell(1,3);

for n = 1:3
Binary{n} = imfill(Binary{n},'holes');
end

for n = 1:3
Dilated{n} = imdilate(Binary{n},strel('disk',15));
end

%%
%Changed active contour replicate value to 300 from 600
%CHanged from 300 to 200
ActiveContour = cell(1,3);
for n = 1:3
ActiveContour{n} = activecontour(Gray{n},Binary{n},300,'Chan-Vese',0.01);
end




%%
%Enlarge areas taken to make sure they cut the outline border of the image
for n = 1:3
   ActiveContour{n} = imdilate(ActiveContour{n},strel('disk',10));
end

%%
%Add missing tissue without color to the equation
for n = 2:3
ActiveContour{n} = ActiveContour{n} + Background_Tissue{n};
ActiveContour{n} = imbinarize(ActiveContour{n});
end

%%
JZ_Removal = (Binary_Outline - ActiveContour{1}) - ActiveContour{3};
%%
%Changed to 10 from 30 to 5 now
JZ_Removal = imdilate(JZ_Removal,strel('disk',5));
JZ_Removal = imfill(JZ_Removal,'holes');

%%
JZ_Removal = imbinarize(JZ_Removal);
%%
[rows, columns, numberOfColorChannels2] = size(JZ_Removal); 
drawnow;
numberOfClusters2 = 2;
grayLevels2 = double(JZ_Removal(:)); 
[clusterIndexes2, clusterCenters2] = kmeans(grayLevels2, numberOfClusters2,...
    'Distance','hamming',...
    'Replicates',1);
labeledImage2 = reshape(clusterIndexes2,rows,columns);
coloredLabels = label2rgb (labeledImage2, 'hsv', 'k', 'shuffle');
[maxValue2, indexOfMaxValue2] = max(clusterCenters2);

JZ_Deduced = labeledImage2 == indexOfMaxValue2;
JZ_Deduced = bwareafilt(JZ_Deduced,1);





%%
%From 15 to 20
JZ_Removal = imdilate(JZ_Deduced,strel('disk',20));
JZ_Removal = Binary_Outline - JZ_Removal;
JZ_Removal = imbinarize(JZ_Removal);

%%
%LZ found by taking largest hamming structure found in image left over after removing JZ
[rows, columns, numberOfColorChannels2] = size(JZ_Removal); 
drawnow;
numberOfClusters2 = 2;
grayLevels2 = double(JZ_Removal(:)); 
[clusterIndexes2, clusterCenters2] = kmeans(grayLevels2, numberOfClusters2,...
    'Distance','hamming',...
    'Replicates',1);
labeledImage2 = reshape(clusterIndexes2,rows,columns);
coloredLabels = label2rgb (labeledImage2, 'hsv', 'k', 'shuffle');
[maxValue2, indexOfMaxValue2] = max(clusterCenters2);

LZ_Deduced = labeledImage2 == indexOfMaxValue2;
LZ_Deduced = bwareafilt(LZ_Deduced,1);




%%
DC_Deduced = (Binary_Outline - LZ_Deduced) - JZ_Deduced;
DC_Deduced = imopen(DC_Deduced,strel('disk',20));

%%
%Outline is the original resized image
Color_LZ = Outline;
Color_JZ = Outline;
Color_DC = Outline;

%%
%Extraction Loop

    LZ_Deduced(:,:,2) = LZ_Deduced;
    LZ_Deduced(:,:,3) = LZ_Deduced(:,:,1);
    Color_LZ(LZ_Deduced == 0) = 0;

    JZ_Deduced(:,:,2) = JZ_Deduced;
    JZ_Deduced(:,:,3) = JZ_Deduced(:,:,1);
    Color_JZ(JZ_Deduced == 0) = 0;

    DC_Deduced(:,:,2) = DC_Deduced;
    DC_Deduced(:,:,3) = DC_Deduced(:,:,1);
    Color_DC(DC_Deduced == 0) = 0;



%%
%End of Zone Isolation
    
    figure,imshow(Color_LZ),title('LZ')
    figure,imshow(Color_JZ),title('JZ')
    figure,imshow(Color_DC),title('DC')
    
    
%%
LZ_Binary = rgb2gray(Color_LZ);    
JZ_Binary = rgb2gray(Color_JZ);
DC_Binary = rgb2gray(Color_DC); 

%Area of Sections
LZ_Binary = imbinarize(LZ_Binary);
JZ_Binary = imbinarize(JZ_Binary);
DC_Binary = imbinarize(DC_Binary);

JZ_Binary_Dilated = imdilate(JZ_Binary,strel('disk',5));
JZ_Background_Perimeter = imfill(JZ_Binary_Dilated,'holes');

LZ_Background = imfill(LZ_Binary,'holes');
JZ_Background_MissingTissue = imfill(JZ_Binary,'holes');
DC_Background = imfill(DC_Binary,'holes');

%Perimeter of Sections
LZ_Perim = bwperim(LZ_Background,8);
JZ_Perim = bwperim(JZ_Background_Perimeter,8);
DC_Perim = bwperim(DC_Background,8);

LZ_Binary_Inv = ~LZ_Binary;
JZ_Binary_Inv = ~JZ_Binary;
DC_Binary_Inv = ~DC_Binary;
LZ_Background = ~LZ_Background;
JZ_Background_MissingTissue = ~JZ_Background_MissingTissue;
DC_Background = ~DC_Background;
LZ_Binary_Inv = LZ_Binary_Inv - LZ_Background; 
JZ_Binary_Inv = JZ_Binary_Inv - JZ_Background_MissingTissue;
DC_Binary_Inv = DC_Binary_Inv - DC_Background;

%Missing tissue in Sections
LZ_Binary_Inv = imbinarize(LZ_Binary_Inv);
JZ_Binary_Inv = imbinarize(JZ_Binary_Inv);
DC_Binary_Inv = imbinarize(DC_Binary_Inv);

%Cross section between LZ and JZ
LZ_Enlarged = imdilate(LZ_Binary,strel('disk',20));
JZ_Enlarged = imdilate(JZ_Binary,strel('disk',20));
Binary_Difference = (double(LZ_Enlarged) - double(JZ_Enlarged)) == 0; 
CrossSection_Background = LZ_Enlarged + JZ_Enlarged;
CrossSection_Background = ~CrossSection_Background;
Binary_Difference = Binary_Difference - CrossSection_Background;

%Cross Section Image
Binary_Difference = imbinarize(Binary_Difference);
Cross_Section_Image = Binary_Difference;

%Value calculated skeleton
Binary_Difference = bwmorph(Binary_Difference,'thin',inf);
%%
%Final Calculations
fprintf('Results \n');

LZ_Size = bwarea(LZ_Binary);
JZ_Size = bwarea(JZ_Binary);
DC_Size = bwarea(DC_Binary);
Voids_LZ = bwarea(LZ_Binary_Inv);
Voids_JZ = bwarea(JZ_Binary_Inv);
Voids_DC = bwarea(DC_Binary_Inv);
PercentMissing_LZ = ((Voids_LZ/(Voids_LZ+LZ_Size))*100);
PercentMissing_JZ = ((Voids_JZ/(Voids_JZ+JZ_Size))*100);
PercentMissing_DC = ((Voids_DC/(Voids_DC+DC_Size))*100);
LZ_Perimeter = bwarea(LZ_Perim);
JZ_Perimeter = bwarea(JZ_Perim);
DC_Perimeter = bwarea(DC_Perim);
LineDistance = bwarea(Binary_Difference);

%%
% Creating New Folder for saving data
mkdir(Title);
BackSlash = '\';
TitleDirectory = [Title BackSlash]; %here
DirectoryName = 'C:\Users\samfi\Desktop\downloads\';        %choose directory for output folder
FileName = [DirectoryName TitleDirectory];


%%
%Writing Data
%Printing Data
CSVExt = '.xlsx';
SpecificName = [Title CSVExt];
FullName = [FileName SpecificName];
%DataName = {'LZ Tissue';'JZ Tissue';'Decidua Tissue';'LZ Missing Tissue';'JZ Missing Tissue';'Decidua Missing Tissue';'LZ Percent Missing';'JZ Percent Missing';'DC Percent Missing';'LZ Perimeter';'JZ Perimeter';'DC Perimeter';'Line Length';};
Data = [LZ_Size;JZ_Size;DC_Size;Voids_LZ;Voids_JZ;Voids_DC;PercentMissing_LZ;PercentMissing_JZ;PercentMissing_DC;LZ_Perimeter;JZ_Perimeter;DC_Perimeter;LineDistance];
%xlswrite(FullName,DataName);
xlswrite(FullName,Data);

%%
%Saving Images

%Segmented Zones
SpecificName = 'LZ Zone Segmentation.jpg';
FullName = [FileName SpecificName];
imwrite(Color_LZ,FullName);

SpecificName = 'JZ Zone Segmentation.jpg';
FullName = [FileName SpecificName];
imwrite(Color_JZ,FullName);

SpecificName = 'Decidua Zone Segmentation.jpg';
FullName = [FileName SpecificName];
imwrite(Color_DC,FullName);

SpecificName = 'Cross Section.jpg';
FullName = [FileName SpecificName];
imwrite(Cross_Section_Image,FullName);

SpecificName = 'Perimeter LZ Region.jpg';
FullName = [FileName SpecificName];
imwrite(LZ_Perim,FullName);

SpecificName = 'Perimeter JZ Region.jpg';
FullName = [FileName SpecificName];
imwrite(JZ_Perim,FullName);

SpecificName = 'Perimeter Decidua Region.jpg';
FullName = [FileName SpecificName];
imwrite(DC_Perim,FullName);

%%
%Finish Stamp
fprintf('Code Ended at %s\n', datestr(now,'HH:MM:SS'))
clc
clearvars -except Image_Range Image_Name File_Name Script_Cycle
end