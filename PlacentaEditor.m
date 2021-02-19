      warning off

File_Question = 'What is the file you are trying to fix? ';
File_Title = input(File_Question,'s');

PlacentaDirectory = 'C:\Users\';
PlacentaFile_Location = [PlacentaDirectory File_Title];
PlacentaSyntax = '\';
Placenta = [PlacentaFile_Location PlacentaSyntax];


fprintf('Which zone are you going to edit \n')


fprintf('1 for Decidua removal from LZ \n')
fprintf('2 for Chorionic Plate removal from JZ \n')
fprintf('3 for Yolk Sack removal from LZ \n')

PlacentaSection_Question = ('What section are you editing?');
PlacentaSection = input(PlacentaSection_Question);

PlacentaJZ = 'JZ Zone Segmentation.jpg';
PlacentaLZ = 'LZ Zone Segmentation.jpg';
PlacentaLZ_Yolk = 'LZ Zone Segmentation.jpg';

if PlacentaSection == 1
    Selected = PlacentaLZ;
    CSV_Subtitle = '_DC_Removal';
elseif PlacentaSection == 2
    Selected = PlacentaJZ;
    CSV_Subtitle = '_Chorionic_Plate_Removal';
elseif PlacentaSection == 3
    Selected = PlacentaLZ_Yolk;
    CSV_Subtitle = '_Yolk_Removal';
end

PlacentaChoice = [Placenta Selected];
SectionEdit = imread(PlacentaChoice);

switch (PlacentaSection)
    
    case 1
JZ_Directory = 'JZ Zone Segmentation.jpg';
PlacentaJZ = [Placenta JZ_Directory];
Placenta_JZ = imread(PlacentaJZ);


figure,imshow(SectionEdit)
h = drawfreehand;
Mask = h.createMask;
close

SectionLZ_Copy = SectionEdit;
SectionDC_Copy = SectionEdit;

Mask_DC = Mask;
Mask_LZ = ~Mask;

Mask_LZ(:,:,2) = Mask_LZ;
Mask_LZ(:,:,3) = Mask_LZ(:,:,1);
SectionLZ_Copy(Mask_LZ == 0) = 0;

Mask_DC(:,:,2) = Mask_DC;
Mask_DC(:,:,3) = Mask_DC(:,:,1);
SectionDC_Copy(Mask_DC == 0) = 0;



CorrectedSection = rgb2gray(SectionLZ_Copy);
CorrectedSection = imbinarize(CorrectedSection);

MaskedSection = rgb2gray(SectionDC_Copy);
MaskedSection = imbinarize(MaskedSection);



LZ_Binary = CorrectedSection;
MaskedSize = bwarea(MaskedSection);


Filled_LZ = imfill(CorrectedSection,'holes');
INV_Filled_LZ = ~Filled_LZ;

VoidSection_LZ = CorrectedSection + INV_Filled_LZ;

LZ_Binary_Inv = ~VoidSection_LZ;

%%

JZ_Binary = rgb2gray(Placenta_JZ);
JZ_Binary = imbinarize(JZ_Binary);
Filled_JZ = imfill(JZ_Binary,'holes');

%Enlarged Areas to create Cross Section Line
Dilated_LZ = imdilate(Filled_LZ,strel('disk',20));
Dilated_JZ = imdilate(Filled_JZ,strel('disk',20));

%Create Background for cross section
JointSections = Dilated_LZ + Dilated_JZ;
JointSections = imfill(JointSections,'holes');
JointSections = ~JointSections;

Binary_Difference = (double(Dilated_LZ) - double(Dilated_JZ)) == 0;
Binary_Difference = Binary_Difference - JointSections;



Binary_Line = bwmorph(Binary_Difference,'thin',inf);

LZ_Perim = bwperim(Filled_LZ,8);
%%
%Calculations
Voids_LZ = bwarea(LZ_Binary_Inv);
LZ_Size = bwarea(LZ_Binary);
LZ_Perimeter = bwarea(LZ_Perim);
PercentMissing_LZ = ((Voids_LZ/(Voids_LZ+LZ_Size))*100);
LineDistance = bwarea(Binary_Line);
Filler = 0;
%%
% Creating New Folder for saving data
mkdir(File_Title);
BackSlash = '\';
TitleDirectory = [File_Title BackSlash]; %here
DirectoryName = 'C:\Users\'; 
FileName = [DirectoryName TitleDirectory];


%%
%Writing Data
%Printing Data
CSVExt = '.xlsx';
SubTitle = [File_Title CSV_Subtitle];
SpecificName = [SubTitle CSVExt];
FullName = [FileName SpecificName];
%DataName = {'LZ Tissue';'JZ Tissue';'Decidua Tissue';'LZ Missing Tissue';'JZ Missing Tissue';'Decidua Missing Tissue';'LZ Percent Missing';'JZ Percent Missing';'DC Percent Missing';'LZ Perimeter';'JZ Perimeter';'DC Perimeter';'Line Length';};
Data = [LZ_Size;Filler;Filler;Voids_LZ;Filler;Filler;PercentMissing_LZ;Filler;Filler;LZ_Perimeter;Filler;Filler;LineDistance];
%xlswrite(FullName,DataName);
xlswrite(FullName,Data);

%%
%Saving Images

%Cross Section Intersection
SpecificName = 'Cross Section Intersection.jpg';
FullName = [FileName SpecificName];
imwrite(Binary_Difference,FullName);

%Segmented Zones
SpecificName = 'LZ Zone Segmentation.jpg';
FullName = [FileName SpecificName];
imwrite(SectionLZ_Copy,FullName);


SpecificName = 'Perimeter LZ Region.jpg';
FullName = [FileName SpecificName];
imwrite(LZ_Perim,FullName);

    case 2
        
LZ_Directory = 'LZ Zone Segmentation.jpg';
PlacentaLZ = [Placenta LZ_Directory];
Placenta_LZ = imread(PlacentaLZ);


figure,imshow(SectionEdit)
h = drawfreehand;
Mask = h.createMask;
close

SectionJZ_Copy = SectionEdit;
SectionCPlate_Copy = SectionEdit;

Mask_CPlate = Mask;
Mask_JZ = ~Mask;

Mask_JZ(:,:,2) = Mask_JZ;
Mask_JZ(:,:,3) = Mask_JZ(:,:,1);
SectionJZ_Copy(Mask_JZ == 0) = 0;

Mask_CPlate(:,:,2) = Mask_CPlate;
Mask_CPlate(:,:,3) = Mask_CPlate(:,:,1);
SectionCPlate_Copy(Mask_CPlate == 0) = 0;


Binary_LZ = rgb2gray(Placenta_LZ);
Binary_LZ = imbinarize(Binary_LZ);
Binary_LZ = imfill(Binary_LZ,'holes');

Binary_JZ = rgb2gray(SectionJZ_Copy);
Binary_JZ = imbinarize(Binary_JZ);
Filled_JZ = imfill(Binary_JZ,'holes');

JZ_Binary_Inv = ~Filled_JZ;

JZ_Binary_Inv = JZ_Binary_Inv + Binary_JZ;
JZ_Binary_Inv = ~JZ_Binary_Inv;

Dilated_JZ_P = imdilate(Binary_JZ,strel('disk',5));
JZ_Background_P = imfill(Dilated_JZ_P,'holes');

Dilated_JZ = imdilate(Filled_JZ,strel('disk',20));
Dilated_LZ = imdilate(Binary_LZ,strel('disk',20));

JointSections = Dilated_JZ + Dilated_LZ;
JointSections = ~JointSections;

Binary_Difference = (double(Dilated_LZ) - double(Dilated_JZ)) == 0;
Binary_Difference = Binary_Difference - JointSections;

Binary_Line = bwmorph(Binary_Difference,'thin',inf);
LZ_Perim = bwperim(Binary_LZ,8);

JZ_Perim = bwperim(JZ_Background_P,8);

%%
%Calculations
JZ_Size = bwarea(Binary_JZ);
Voids_JZ = bwarea(JZ_Binary_Inv);
PercentMissing_JZ = ((Voids_JZ/(Voids_JZ+JZ_Size))*100);
JZ_Perimeter = bwarea(JZ_Perim);
LineDistance = bwarea(Binary_Line);
Filler = 0;
%%
% Creating New Folder for saving data
mkdir(File_Title);
BackSlash = '\';
TitleDirectory = [File_Title BackSlash]; %here
DirectoryName = 'C:\Users\'; 
FileName = [DirectoryName TitleDirectory];


%%
%Writing Data
%Printing Data
CSVExt = '.xlsx';
SubTitle = [File_Title CSV_Subtitle];
SpecificName = [SubTitle CSVExt];
FullName = [FileName SpecificName];
%DataName = {'LZ Tissue';'JZ Tissue';'Decidua Tissue';'LZ Missing Tissue';'JZ Missing Tissue';'Decidua Missing Tissue';'LZ Percent Missing';'JZ Percent Missing';'DC Percent Missing';'LZ Perimeter';'JZ Perimeter';'DC Perimeter';'Line Length'};
Data = [Filler;JZ_Size;Filler;Filler;Voids_JZ;Filler;Filler;PercentMissing_JZ;Filler;Filler;JZ_Perimeter;Filler;LineDistance];
%xlswrite(FullName,DataName);
xlswrite(FullName,Data);

%%
%Saving Images

%Cross Section Intersection
SpecificName = 'Cross Section Intersection.jpg';
FullName = [FileName SpecificName];
imwrite(Binary_Difference,FullName);

%Segmented Zones
SpecificName = 'JZ Zone Segmentation.jpg';
FullName = [FileName SpecificName];
imwrite(SectionJZ_Copy,FullName);


SpecificName = 'Perimeter JZ Region.jpg';
FullName = [FileName SpecificName];
imwrite(JZ_Perim,FullName);


    case 3
figure,imshow(SectionEdit)
h = drawfreehand;
Mask = h.createMask;
close

SectionLZ_Copy = SectionEdit;

Mask_LZ = ~Mask;

Mask_LZ(:,:,2) = Mask_LZ;
Mask_LZ(:,:,3) = Mask_LZ(:,:,1);
SectionLZ_Copy(Mask_LZ == 0) = 0;

Binary_LZ = rgb2gray(SectionLZ_Copy);
Binary_LZ  = imbinarize(Binary_LZ);

Filled_LZ = imfill(Binary_LZ,'holes');

Binary_LZ_Inv = ~Filled_LZ;

Binary_LZ_Inv = Binary_LZ_Inv + Binary_LZ;
Binary_LZ_Inv = ~Binary_LZ_Inv;

LZ_Perim = bwperim(Filled_LZ,8);
%%
%Calculations
Voids_LZ = bwarea(Binary_LZ_Inv);
LZ_Size = bwarea(Binary_LZ);
LZ_Perimeter = bwarea(LZ_Perim);
PercentMissing_LZ = ((Voids_LZ/(Voids_LZ+LZ_Size))*100);
Filler = 0;

%%
%Printing

%Creating New Folder for saving data
mkdir(File_Title);
BackSlash = '\';
TitleDirectory = [File_Title BackSlash]; %here
DirectoryName = 'C:\Users\'; 
FileName = [DirectoryName TitleDirectory];


%%
%Writing Data
%Printing Data
CSVExt = '.xlsx';
SubTitle = [File_Title CSV_Subtitle];
SpecificName = [SubTitle CSVExt];
FullName = [FileName SpecificName];
%DataName = {'LZ Tissue';'JZ Tissue';'Decidua Tissue';'LZ Missing Tissue';'JZ Missing Tissue';'Decidua Missing Tissue';'LZ Percent Missing';'JZ Percent Missing';'DC Percent Missing';'LZ Perimeter';'JZ Perimeter';'DC Perimeter';'Line Length';};
Data = [LZ_Size;Filler;Filler;Voids_LZ;Filler;Filler;PercentMissing_LZ;Filler;Filler;LZ_Perimeter;Filler;Filler;Filler];
%xlswrite(FullName,DataName);
xlswrite(FullName,Data);

%%

%Segmented Zones
SpecificName = 'LZ Zone Segmentation.jpg';
FullName = [FileName SpecificName];
imwrite(SectionLZ_Copy,FullName);

SpecificName = 'Perimeter LZ Region.jpg';
FullName = [FileName SpecificName];
imwrite(LZ_Perim,FullName);
end


clearvars
clc
