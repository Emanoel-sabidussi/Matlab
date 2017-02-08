clear all
clc
close all

x = zeros(50,8);
BaseName='L';
extName = '.csv';
t = (linspace(1,50,50))';
for k=1:50
    FileName = [BaseName,num2str(k)];
    FileName = strcat(FileName, extName);
   
    x(k,:) = csvread(FileName,0,0);
    xS(k,:) = sqrt(x(k,:));
    xS(k,:);
    
   
end

z = zeros(50,8);
BaseName='';
extName = '.csv';
t = (linspace(1,50,50))';
figure
for k=1:50
    FileName = [BaseName,num2str(k)];
    FileName = strcat(FileName, extName);
   
    z(k,:) = csvread(FileName,0,0);
    zS(k,:) = sqrt(z(k,:));
    zS(k,:);
    
   
end

%%

for k=1:8
    scatter (t, xS(:,k));
    hold on
    
end

legend('Mean','STD','Energy','Entropy','ASM','EntropyGLCM','Correlation','IDM')
title ('Correlation x IDM');


X = [xS; zS];
Y = zeros(50,1);
Y1 = ones(50,1);
Y = [Y ; Y1];


%%
%clear all
close all

pca_t = linspace(1,50,50);
pca_t1 = linspace(1,30,30);
PCA_Proj_Result = csvread('PCA_Projection_Result.csv');
PCA_Eigen_Vectors = csvread('PCA_Eigen_Vectors.csv');
PCA_Back_Project = csvread('PCA_Back_Project.csv');
Feature_Vector = csvread('Feature_Vector.csv');


for k =1:8
    scatter(pca_t,Feature_Vector(k,1:50),'bo');
    hold on
end

for k =1:8
    scatter(pca_t,Feature_Vector(k,51:100),'ro');
    hold on
end

for k =1:8
    scatter(pca_t1,Feature_Vector(k,101:130),'go');
    hold on
end

legend('Layer1','Layer2','Layer3');
title('Feature_Vector');



for i=1:7
    for j=2:8
   if i~=j
        figure
        scatter(Feature_Vector(i,1:50),Feature_Vector(j,1:50),'bo');
        hold on
        scatter(Feature_Vector(i,51:100),Feature_Vector(j,51:100),'ro');
        hold on
        scatter(Feature_Vector(i,101:130),Feature_Vector(j,101:130),'go');
        legend('Layer1','Layer2','Layer3');
        str1 = num2str(i);
        str2 = num2str(j);
        str2 = strcat(str1, str2);
        title(str2);
   end
    end
end

figure
scatter(Feature_Vector(2,1:50),Feature_Vector(4,1:50),'bo');
hold on
scatter(Feature_Vector(2,51:100),Feature_Vector(4,51:100),'ro');
hold on
scatter(Feature_Vector(2,101:130),Feature_Vector(4,101:130),'go');
legend('Layer1','Layer2','Layer3');
title('STD x Energy');



figure
for k = 1:4
    scatter(pca_t,PCA_Proj_Result(k,1:50),'bo');
    hold on
end

for k = 1:4
    scatter(pca_t,PCA_Proj_Result(k,51:100),'ro');
    hold on
end

for k = 1:4
    scatter(pca_t1,PCA_Proj_Result(k,101:130),'go');
    hold on
end

title('PROJECT_RESULT');

figure
scatter(PCA_Proj_Result(1,1:50),PCA_Proj_Result(2,1:50),'bo');
hold on
scatter(PCA_Proj_Result(1,51:100),PCA_Proj_Result(2,51:100),'ro');
hold on
scatter(PCA_Proj_Result(1,101:130),PCA_Proj_Result(2,101:130),'go');
legend('Layer1','Layer2');
title('Feat2 x Feat3');





pca_t = linspace(1,9,9);

figure
for k = 1:4
   
    scatter(pca_t,PCA_Eigen_Vectors(k,:));
    hold on
    
end
title('PCA_Eigen_Vectors');

pca_t = linspace(1,50,50);
figure
for k = 1:9
   
    scatter(pca_t,PCA_Back_Project(k,1:50),'bo');
    hold on
    
end

for k = 1:9
   
    scatter(pca_t,PCA_Back_Project(k,51:100),'ro');
    hold on
    
end

for k = 1:9
   
    scatter(pca_t1,PCA_Back_Project(k,101:130),'go');
    hold on
    
end
title('PCA_Back_Project');




%%

figure
scatter3(Feature_Vector(2,1:50),Feature_Vector(4,1:50),Feature_Vector(8,1:50),'bo');
hold on
scatter3(Feature_Vector(2,51:100),Feature_Vector(4,51:100),Feature_Vector(8,51:100),'ro');

%%



SVMModel = fitcsvm(X,Y,'BoxConstraint',100, 'Standardize', true,'DeltaGradientTolerance',1e-1, 'KernelFunction','RBF','KernelScale','auto');

classOrder = SVMModel.ClassNames
sv = SVMModel.SupportVectors;
figure
scatter3(xS(:,4), xS(:,7),xS(:,2))
hold on
scatter3(zS(:,4), zS(:,7), zS(:,2))
hold on
plot3(sv(:,4),sv(:,7),sv(:,2), 'ko', 'MarkerSize', 10)










%%
legend('Mean','STD','Energy','Entropy','ASM','EntropyGLCM','Correlation','IDM')
xlabel('Samples') % x-axis label
ylabel('Features') % y-axis label
%%
for k=1:7
    %figure
    %scatter(xS(k,:),xS(k+1,:));
    
    titl=[BaseName,num2str(k)];
    titl1=[BaseName,num2str(k+1)];
    
    titl = [titl titl1];
    title(titl);
end










%%
l = zeros(8,5);
figure
BaseName='l';
extName = '.csv';
t = linspace(1,8,8)';
for k=1:5
    
    FileName = strcat(FileName, extName);
   
    l(:,k) = csvread(FileName,0,0);
    lS(:,k) = sqrt(l(:,k));
    
    scatter(lS(5,k),lS(6,k));
    hold on
    
end

GLCM = csvread('GLCM.csv', 0, 0 );
GLCM = GLCM.^2;

result = sum(sum(GLCM));




