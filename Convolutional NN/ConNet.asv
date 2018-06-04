clear %clear any variables from previous runs
clc
%29.9408
filtersize=3; %size of each convolutional filter, must be odd
dx=.01; %step size used in central distance formula for calculating the loss gradient during training
precision = .1 %desired precision of network outputs
step = .000000005/(filtersize*filtersize); %step size for changing the filter value

path=pwd;
trueTemp = load('tempvals');  %load in the temp matrix
data = load('imagevals');

[height,length]=size(data);

                                %this section adds a layer of zeros around the border of the data in order
                                %to allow the filters to convolute over
                                %every location in the data.  future coding
                                %will use height, length, etc as if the
                                %padding did not exist, but the padsize
                                %variable will be added to indexes as
padsize = (filtersize-1)/2;     %needed to fetch and set correct locations.
newMat = zeros(height+2*padsize,length+2*padsize);
for m=1:height
    for n=1:length
        newMat(m+padsize,n+padsize)=data(m,n);
    end
end
data = newMat;
[newHeight,newLength]=size(data);

%conFilterNum = ceil((sqrt(height^2 + length^2))/2); comments are because
%we are using the filter-sharing assumption
conFilter1=zeros(filtersize,filtersize...,conFilterNum... 
    );
%for n1=1:conFilterNum
n = height*length;
    for n1=1:filtersize
        for n2=1:filtersize
            conFilter1(n1,n2) = randn() / sqrt(n);
        end
    end
%end
%create matrix of biases as well?


loss = 0;
loss2 = zeros(filtersize,filtersize);
loss3 = zeros(filtersize,filtersize);
grad = zeros(filtersize,filtersize);
sum1 = 0;
conresults = zeros(height,length);
for x=1:300 %number of times to run for training
    
    loss = 0;
    loss2 = zeros(filtersize,filtersize);
    loss3 = zeros(filtersize,filtersize);
    
    vals1=zeros(filtersize,filtersize); %values are initially unsummed in order to calculate loss for each
    vals2 = zeros(filtersize,filtersize); %uses +dx for Central Difference Theorem
    vals3 = zeros(filtersize,filtersize); %uses -dx for CDT
    sum1=0;
    for m=1+padsize:height-padsize
        for n=1+padsize:length-padsize
            vals1=zeros(filtersize,filtersize);
            vals2 = zeros(filtersize,filtersize); 
            vals3 = zeros(filtersize,filtersize);
            for f1=1:filtersize
                for f2=1:filtersize
                    vals1(f1,f2) =(conFilter1(f1,f2,1 ...CHANGE IF GETTING RID OF FILTER SHARING ...
                        ) * data(m+f1-1,n+f2-1)); %compute summed dot product
                    vals2(f1,f2) =((conFilter1(f1,f2,1) + dx) * data(m+f1-1,n+f2-1));  
                    vals3(f1,f2) =((conFilter1(f1,f2,1) - dx) * data(m+f1-1,n+f2-1)); %of filter and filtersize x
                                                                                    %filtersize region centered
                                                                                    %about data(mxn)
                end
            end
            
            sum1 = sum(sum(vals1));  %total loss is determined for the purpose of figuring out whether
            conresults(m,n) = sum1;                           %to keep training, other losses are
                                       %for determining gradient of each
                                       %filter value
            loss = loss + (sum1 - trueTemp(m,n))^2;
            for f1=1:filtersize
                for f2=1:filtersize
                    loss2(f1,f2) = loss2(f1,f2) + ((sum1 - vals1(f1,f2) + vals2(f1,f2))- trueTemp(m,n))^2;
                    loss3(f1,f2) = loss3(f1,f2) + ((sum1 - vals1(f1,f2) + vals3(f1,f2))- trueTemp(m,n))^2;
                end
            end
            
            %sum = 0; <-- do it this way for application if not training
            %for f1=1:filtersize
            %    for f2=1:filtersize
            %        sum = sum + ((conFilter1(f1,f2,1 ...CHANGE IF GETTING RID OF FILTER SHARING ...
            %            )) * data(m+f1-1,n+f2-1)); %compute summed dot product
            %                                                                        %of filter and filtersize x
            %                                                                        %filtersize region centered
            %                                                                        %about data(mxn)
            %    end
            %end
        end
    end
    %calculate loss
 
    loss = loss/(height*length)
    %check if close enough
    if loss<precision
        break
    end
    %change filters based on gradient
    loss2 = loss2./(height*length);
    loss3 = loss3./(height*length);
    grad = loss2 - loss3;
    conFilter1 = conFilter1 - (grad .* step);
    
    
end
%REMOVED SECTION: we are using a filter sharing assumption instead of a
%concentric filter system.
%centerm = height/2;
%centern = length/2;
%dist = 0;
%mdist = 0;
%ndist = 0;
%filterNum = 0;
%filterVectors = zeros(conFilterNum-1,conFilterNum*4);
%for n1=1:conFilterNum  %because matlab doesn't do well as storing a matrix of diffenently sized matrices, 
                        %the first element for each filter vector will list
                        %the index for where to put the next element.
                        %Similarly, the size of each vector may be computed
                        %as filterVectors(n1,1) - 2.
%    filterVectors(n1,1)=2;
%end
%index2 = 0;
%sum = 0;
%for m=1:height
%    m
%    for n=1:length
        %mdist = centerm - m + .5;
        %ndist = centern - n + .5;
        %dist = sqrt(mdist^2 + ndist^2);
%        filterNum = floor(dist) + 1;
%        index2 = filterVectors(filterNum,1);
%        sum = 0;
%        for f1=1:filtersize
%            for f2=1:filtersize
%                sum = sum + (conFilter(f1,f2,filterNum) * data(m+f1-2+padsize,n+f2-2+padsize)); %compute summed dot product
                                                                                %of filter and filtersize x
                                                                                %filtersize region centered
                                                                                %about data(mxn)
%            end
%        end
%        filterVectors(filterNum,index2)=sum;
%        filterVectors(filterNum,1)= filterVectors(filterNum,1)+1;
%    end
%end