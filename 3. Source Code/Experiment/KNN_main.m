clear;
close all;
tic;
%% Variables

load('Image_labels.mat');
load('Boundaries_of_check_region.mat');
image_row=140;
image_col=140;
region_x=10;
region_y=10;
addpath(genpath('./Feature_Extraction/'));
image_path='./Cropped_images_JAFFE/';
img_list=dir([image_path,'*.tiff']);
Feature_Histograms=zeros(size(img_list,1),region_x*region_y*(256+256+2+88)); %proposed

%% Feature Vector Generation
disp('Feature generating...');
load('Feature_Histograms.mat');
for i=1:size(img_list,1)
    fprintf('Accessing Image# %d..\n',i);
    img=imresize(imread([image_path,img_list(i).name]),[image_row image_col ]);   
    left_eyebrow=Boundaries_of_check_region(i,1);
    right_eyebrow=Boundaries_of_check_region(i,2);
    lower_eye=Boundaries_of_check_region(i,3);
    upper_lip=Boundaries_of_check_region(i,4);    
    Feature_Histograms(i,:)=proposed_method_exp2(img,left_eyebrow,right_eyebrow,lower_eye,upper_lip, region_x,region_y);
    
end
save('Feature_Histograms.mat','Feature_Histograms');

%% Classification
Num_of_fold=5;
ValueofK_for_KNN=1;
indices = crossvalind('Kfold',Image_labels,10);

sum_accuracy=0;
min_accuracy=1;
max_accuracy=0;
for i = 1:Num_of_fold
    test_indices = (indices == i); train_indices = ~test_indices;  
    %% KNN
    disp('KNN Training...');
    knn = fitcknn(Feature_Histograms(train_indices,:), Image_labels(train_indices),'NumNeighbors',ValueofK_for_KNN);
    disp('KNN Testing...');
    experiment_result = predict(knn, Feature_Histograms(test_indices,:));
      
    %% Performance
    cp = classperf(Image_labels);
    classperf(cp,experiment_result,test_indices);
    fprintf('Iteration # %d =>',i);
    
    if (cp.CorrectRate<min_accuracy)
        min_accuracy=cp.CorrectRate;
    end
    if (cp.CorrectRate>max_accuracy)
        max_accuracy=cp.CorrectRate;
    end    
    
    sum_accuracy=sum_accuracy+cp.CorrectRate;
    fprintf('Accuracy: %f \n',cp.CorrectRate*100);
    
    %% Confusion Matrix
    disp('Confusion Matrix');
    C=confusionmat(Image_labels(test_indices),experiment_result);
    disp(C);

    load('Unique_Class.mat');
    Num_of_class=size(Unique_Class,1);
    Num_of_Test_img=sum(test_indices);
    Targets=zeros(Num_of_class,Num_of_Test_img);
    Outputs=zeros(Num_of_class,Num_of_Test_img);
    for m=1:Num_of_class
        Targets(m,:)=(Image_labels(test_indices)==Unique_Class(m));
        Outputs(m,:)=(experiment_result==Unique_Class(m));
    end
    figure;
    plotconfusion(Targets,Outputs,sprintf('Iteration #%d ',i));

end
avg_accuracy=sum_accuracy/Num_of_fold;

fprintf('\nMax Accuracy: %f \n',max_accuracy*100);
fprintf('min Accuracy: %f \n',min_accuracy*100);
fprintf('Average Accuracy: %f \n',avg_accuracy*100);


%% Write to  Result.txt File
fileID=fopen('Result.txt','at');
time_and_date=datetime;
long_line='---------------------------';
Dataset='JAFFE';
Method='Proposed Method (median_edge_res)';
Additional_info='5 fold'; 
%For KNN
fprintf(fileID,'\n%s\n%s\nKNN\n%s\nValue of K for KNN: %d\nMethod: %s\nDataset: %s\nRegion: %dx%d\nMax Axxuracy: %f\nMin Accuracy: %f\nAvg Accuracy: %f\n%s\n\n\n',time_and_date,long_line,long_line,ValueofK_for_KNN,Method,Dataset,region_x,region_y,max_accuracy*100,min_accuracy*100,avg_accuracy*100,Additional_info);

fclose('all');
toc;

%% Notification sound
load chirp               % handel,gong,laughter,train ,splat
sound(y,Fs)


