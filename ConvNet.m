path=pwd;
%load('grey');
%loading in the images and their associated temperature values
load('temps');
load('rvals');
load('gvals');
load('bvals');
imagemax = 10;
I = zeros(imagemax,1080,1920,3); %storing images in a single matrix
I(:,:,:,1) = Ired;
I(:,:,:,2) = Igreen;
I(:,:,:,3) = Iblue;
c=15; %size of smaller images that we want to break the larger ones into
images= zeros(c,c,3,imagemax*5184*20*20/(c*c));
i=1;
index = randperm(imagemax*5184*20*20/(c*c)); %creating a random index so that images
temps = zeros(1);                               %aren't in a particular order
for n = 1:imagemax %loop breaks larger images into smaller, associates each with 
for x=1:54*(20/c)       %an average temperature
   for y=1:96*(20/c)
      images(:,:,:,index(i)) = I(n,(x*c-(c-1)):(x*c),(y*c-(c-1)):(y*c),:); 
      temps(index(i)) = mean(mean(True_Temp(n,(x*c-(c-1)):(x*c),(y*c-(c-1)):(y*c))));
      i = i+1;
   end
end
end
temps = temps.';
%load('greyval');
load('redtest'); %repeat the process with another dataset for validation
load('greentest');
load('bluetest');
load('tempsval');
I(:,:,:,1) = Ired;
I(:,:,:,2) = Igreen;
I(:,:,:,3) = Iblue;
valimages= zeros(c,c,3,5184*20*20/(c*c));
i = 1;
valtemps = zeros(1);
for x=1:54*(20/c)
   for y=1:96*(20/c)
      n =randi(imagemax);
      valimages(:,:,:,i) = I(n,(x*c-(c-1)):(x*c),(y*c-(c-1)):(y*c),:);
      valtemps(i) = mean(mean(True_Temp(n,(x*c-(c-1)):(x*c),(y*c-(c-1)):(y*c))));
      i = i+1;
   end
end

valtemps = valtemps.';
layers = [ %convolutional layer structure
    imageInputLayer([c c 3])

    convolution2dLayer(3,8,'Padding','same')
    batchNormalizationLayer
    reluLayer
    averagePooling2dLayer(2,'Stride',2)

    convolution2dLayer(5,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    convolution2dLayer(5,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    averagePooling2dLayer(2,'Stride',2)

    convolution2dLayer(6,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
        
    convolution2dLayer(6,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    convolution2dLayer(6,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
   
    dropoutLayer(0.1)
    fullyConnectedLayer(1)
    regressionLayer];

%various options and settings for the CNN.  Many of these are from a
%tutorial on the matlab website. 
validationFrequency = floor(numel(temps)/miniBatchSize);
options = trainingOptions('sgdm',...
    'MiniBatchSize',256,...
    'MaxEpochs',30,...
    'InitialLearnRate',1e-4,...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.1,...
    'LearnRateDropPeriod',20,...
    'Shuffle','every-epoch',...
    'ValidationData',{valimages,valtemps},...
    'ValidationFrequency',validationFrequency,...
    'ValidationPatience',Inf,...
    'Plots','training-progress',...
    'Verbose',false);

%train the network
net = trainNetwork(images,temps,layers,options);