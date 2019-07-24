path=pwd;
%load('grey');
%loading in the images and their associated temperature values
load('TempList');
load('RedList');
load('GreenList');
load('BlueList');
imagemax = 10;
I = zeros(imagemax,1080,1920,3); %storing images in a single matrix
I(:,:,:,1) = Ired; 
I(:,:,:,2) = Igreen;
I(:,:,:,3) = Iblue;
c=60; %size of smaller images that we want to break the larger ones into
images= zeros(c,c,3,imagemax*5184*20*20/(c*c));
i=1;
index = randperm(imagemax*5184*20*20/(c*c)); %creating a random index so that images
temps = zeros(c*c,54*96*(20/c)*(20/c));           %aren't in a particular order
matrix= zeros(c,c,54*96*(20/c)*(20/c));
for n = 1:imagemax %loop breaks larger images into smaller, associates each with 
for x=1:54*(20/c)       %an set of temperatures
   for y=1:96*(20/c)
      images(:,:,:,index(i)) = I(n,(x*c-(c-1)):(x*c),(y*c-(c-1)):(y*c),:); 
      matrix=True_Temp(n,(x*c-(c-1)):(x*c),(y*c-(c-1)):(y*c));
      temps(:,index(i)) = matrix(:);
      i = i+1;
   end
end
end
temps = temps.';
%load('greyval');
load('ttest'); %repeat the process with another dataset for validation
load('btest');
load('gtest');
load('rtest');
I(:,:,:,1) = Ired;
I(:,:,:,2) = Igreen;
I(:,:,:,3) = Iblue;
valimages= zeros(c,c,3,5184*20*20/(c*c));
i = 1;
valtemps = zeros(c*c,54*96*(20/c)*(20/c));
valmatrix=zeros(c,c,54*96*(20/c)*(20/c));
for x=1:54*(20/c)
   for y=1:96*(20/c)
      n =randi(imagemax);
      valimages(:,:,:,i) = I(n,(x*c-(c-1)):(x*c),(y*c-(c-1)):(y*c),:);
      valmatrix = True_Temp(n,(x*c-(c-1)):(x*c),(y*c-(c-1)):(y*c));
      valtemps(:,i) = valmatrix(:);
      i = i+1;
   end
end

valtemps = valtemps.';
 layer1 = [    
     imageInputLayer([c,c,3], 'Name', 'input')
 
     convolution2dLayer(5,8,'Padding','same','Name', 'conv_1')
     batchNormalizationLayer('Name', 'b_1')
     reluLayer('Name','relu_1')
     averagePooling2dLayer(2,'Stride',2,'Name', 'avp_1')
 
     convolution2dLayer(7,64,'Padding','same','Name','conv_2')
     batchNormalizationLayer('Name', 'b_2')
     reluLayer('Name','relu_2')
     
     convolution2dLayer(7,16,'Padding','same','Name','conv_3')
     batchNormalizationLayer('Name','b_3')
     reluLayer('Name','relu_3')
     averagePooling2dLayer(2,'Stride',2,'Name','avp_2')
     
%      convolution2dLayer(7,16,'Padding','same','Name','conv_3x')
%      batchNormalizationLayer('Name','b_3x')
%      reluLayer('Name','relu_3x')
%      averagePooling2dLayer(2,'Stride',2,'Name','avp_2x')
 
     convolution2dLayer(7,32,'Padding','same','Name','conv_4')
     batchNormalizationLayer('Name','b_4')
     reluLayer('Name','relu_4')
         
     convolution2dLayer(6,32,'Padding','same','Name','conv_5')
     batchNormalizationLayer('Name','b_5')
     reluLayer('Name','relu_5')
     
     convolution2dLayer(6,32,'Padding','same','Name','conv_6')
     batchNormalizationLayer('Name','b_6')
     reluLayer('Name','relu_6')
     transposedConv2dLayer(1,32,'Name','tc_1')
    ];
 
% % This is the collection of layers after the first concatenation and before
% % the second concatenation.
 layer2 = [
     convolution2dLayer(5,16,'Padding','same','Name','conv_7')
     batchNormalizationLayer('Name','b_7')
     reluLayer('Name','relu_7')
%      convolution2dLayer(5,16,'Padding','same','Name','conv_x')
%      batchNormalizationLayer('Name','b_x')
%      reluLayer('Name','relu_x')
%      convolution2dLayer(5,16,'Padding','same','Name','conv_y')
%      batchNormalizationLayer('Name','b_y')
%      reluLayer('Name','relu_y')
     transposedConv2dLayer(16,64,'Name','tc_2')
     
 ];

%  layer2b = [
%      convolution2dLayer(5,16,'Padding','same','Name','conv_7b')
%      batchNormalizationLayer('Name','b_7b')
%      reluLayer('Name','relu_7b')
% %      convolution2dLayer(5,16,'Padding','same','Name','conv_x')
% %      batchNormalizationLayer('Name','b_x')
% %      reluLayer('Name','relu_x')
% %      convolution2dLayer(5,16,'Padding','same','Name','conv_y')
% %      batchNormalizationLayer('Name','b_y')
% %      reluLayer('Name','relu_y')
%      transposedConv2dLayer(1,64,'Name','tc_2b')
%      
%  ];

%This is the collection of layers after the second concatenation
layer3 = [
    %imageInputLayer([c,c,3], 'Name', 'input')
    convolution2dLayer(3,16,'Padding','same','Name','conv_8')
    batchNormalizationLayer('Name','b_8')
    reluLayer('Name','relu_8')
%     convolution2dLayer(3,16,'Padding','same','Name','conv_8x')
%     batchNormalizationLayer('Name','b_8x')
%     reluLayer('Name','relu_8x')
%     convolution2dLayer(3,16,'Padding','same','Name','conv_8y')
%     batchNormalizationLayer('Name','b_8y')
%     reluLayer('Name','relu_8y')
%   convolution2dLayer(3,16,'Padding','same','Name','conv_9')
 %   batchNormalizationLayer('Name','b_9')
  %  reluLayer('Name','relu_9')
 %   convolution2dLayer(3,16,'Padding','same','Name','conv_10')
  %  batchNormalizationLayer('Name','b_10')
   % reluLayer('Name','relu_10')
   % convolution2dLayer(3,16,'Padding','same','Name','conv_11')
   % batchNormalizationLayer('Name','b_11')
   % reluLayer('Name','relu_11')
    dropoutLayer(0.1,'Name','drop_1')
    fullyConnectedLayer(c*c,'Name','fullcon_1')
    regressionLayer('Name','reglay_1')
    
];

% To concatenate the two layers you must use the a concat layer
% The first number is the number of outputs it is concating together
% In our U-net this will be 2. 
concat1 = depthConcatenationLayer(2,'Name','concat_1');
concat2 = depthConcatenationLayer(2,'Name','concat_2');
%concat3 = depthConcatenationLayer(2,'Name','concat_3');



% To create a layer graph you must first initialize it like so:
lgraph = layerGraph;

% To add the layers you use the add Layers command
% First input is the layergraph
% Second input is the layer
lgraph = addLayers(lgraph, layer1);
lgraph = addLayers(lgraph, layer2);
%lgraph = addLayers(lgraph, layer2b);
lgraph = addLayers(lgraph, layer3);


%To add a single layer instead of a group use the layer.layer name 
% not the layer name. ie use concat1 not 'concat_1'
lgraph = addLayers(lgraph, concat1);
lgraph = addLayers(lgraph, concat2);
%lgraph = addLayers(lgraph, concat3);

%Here we are connecting the layers.
%With connectLayers command first input is the layer graph
%Second input is the starting layer
%Third input layer that the connection is going to
% with the concat/in1 we are telling the computer that this is the 
% first of the 2 inputs we set it up to have
lgraph = connectLayers(lgraph, 'tc_1', 'concat_1/in1');
lgraph = connectLayers(lgraph, 'relu_4', 'concat_1/in2');
lgraph = connectLayers(lgraph, 'concat_1','conv_7');

lgraph = connectLayers(lgraph, 'tc_2', 'concat_2/in1');
lgraph = connectLayers(lgraph, 'relu_2', 'concat_2/in2');
lgraph = connectLayers(lgraph, 'concat_2', 'conv_8');

% lgraph = connectLayers(lgraph, 'tc_2b', 'concat_3/in1');
% lgraph = connectLayers(lgraph, 'relu_3', 'concat_3/in2');
% lgraph = connectLayers(lgraph, 'concat_3', 'conv_8');
%various options and settings for the CNN.  Many of these are from a
%tutorial on the matlab website. 

miniBatchSize = 256;
validationFrequency = floor(numel(temps)/miniBatchSize);
options = trainingOptions('sgdm',...
    'MiniBatchSize',512,...
    'MaxEpochs',90,...
    'InitialLearnRate',.2,...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',.7,...
    'LearnRateDropPeriod',10,...
    'Shuffle','every-epoch',...
    'ValidationData',{valimages,valtemps},...
    'ValidationFrequency',validationFrequency,...
    'ValidationPatience',Inf,...
    'Plots','training-progress',...
    'Verbose',false);
    
temps = temps;
%train the network
net = trainNetwork(images,temps,lgraph,options);
