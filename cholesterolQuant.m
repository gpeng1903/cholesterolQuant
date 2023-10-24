%% cholesterol count in cultured cells, count GFP/filipin, and recruit intensity etc
% auto read all confocal files in the folder
% import raw confocal data with bio-formats tools as dataRaw
% from dataRaw find sharp/focused slices 
% adaptthresh, generate mask, and count GFP/filipin intensity, and recruit intensity etc
% output result as Excel
% Gang Peng @ 2021


%% pre-processing
% exp. specific setup
thresh_sensitivity = 0.4; % empirically determined for each batch of data, usually 0.3-0.4.
dataChannels = 2; % GFP and filipin

% select slices to count
focusedN = 4; % number of the sharpest/focused slices in each stack = focusedN

% data folder
dataPath = pwd;
dataList = dir(strcat(dataPath, '\*.oib'));
dataN = length(dataList);


% calculation tables
resultFileName = 'countResult.xlsx';



%% calculation loop

dataCountTable = zeros(dataN, 30);

for iData = 1:dataN
    
    % read image stack
    dataName = dataList(iData).name;
    dataRaw = bfopen(dataName); % use bfopen from https://www.openmicroscopy.org/bio-formats/downloads/    
    % dataRaw is a composite cell, containing raw data and annotatons
    
    dataImage = dataRaw{1,1}(:,1); % dataImage is the image stack, in XY-C-Z order
    sliceN = size(dataImage, 1); % get numer of slices in the stack
    
    % check and select confocal slices that are in focus and sharp
    focusT_blue = zeros(sliceN/dataChannels, 1); % foscusT store the focusMeasure value
    sliceI = 1;
    for i = 1:dataChannels:sliceN
        a1 = dataImage{i, 1};
        FM = fmeasure(a1, 'TENV'); % use the TENV method for microscopical images
        focusT_blue(sliceI, 1) = FM;
        % disp(sliceI);
        sliceI = sliceI + 1;    
    end
    % figure, plot(focusT_green, '-o');
    % sort
    [dummy1, idxBlue] = sort(focusT_blue, 'descend');    

    % focus index for green/GFP
    focusT_green = zeros(sliceN/dataChannels, 1);
    sliceI = 1;
    for i = 2:dataChannels:sliceN
        a1 = dataImage{i, 1};
        FM = fmeasure(a1, 'TENV');
        focusT_green(sliceI, 1) = FM;
        % disp(sliceI);
        sliceI = sliceI + 1;    
    end
    % figure, plot(focusT_green, '-o');
    [dummy, idxGreen] = sort(focusT_green, 'descend');
    
    
    
    % count filipin total, background correction then count
    % note filipin is the 1,3,5th... slices          
    % initiate temp array    
    pBlueCount = zeros(1, focusedN); % filipin count
    pBlueArea = zeros(1, focusedN); % filipin area    
    pBlueLevel = zeros(1, focusedN); % slice level determined by adaptthresh 
    
    for sP = 1:focusedN
        % use idxBlue find which slice to use
        bluePickNum = idxBlue(sP)*dataChannels - 1;
        pickedBlueSlice = dataImage{bluePickNum, 1}; % choose focused slice        
        % level and segment
        levelBlue = adaptthresh(pickedBlueSlice, thresh_sensitivity);  
        bw_blue = imbinarize(pickedBlueSlice, levelBlue);
%         figure, imshow(bw);
        bw2_blue = bwareaopen(bw_blue, 3); % remove random noise 
%         figure, imshow(bw2);
       % remove background, use uniform averaged level
       blueBkg = mean(levelBlue(:)) * (2^16); % convert level to uint16 value
       pickBlueRvBkg = imsubtract(pickedBlueSlice, blueBkg); 
       
       % calculate results
        pBlueCount(1, sP) = sum(pickBlueRvBkg(bw2_blue)); % count the blue channel total
        pBlueArea(1, sP) = sum(bw2_blue(:)); % count blue mask area
        pBlueLevel(1, sP) = mean(bw2_blue(:));    
    end    
    dataCountTable(iData, 37:40) = pBlueCount;
    dataCountTable(iData, 43:46) = pBlueArea;
    dataCountTable(iData, 49:52) = pBlueLevel;
    
    % major calculataion below
    % count GFP positive total, use formula:  blue(gfp+ position) / (area of gfp+)
    % also count blue total
    % 
    pBlueOnGreen = zeros(1, focusedN); %
    pBlueTotal = zeros(1, focusedN); %
    pGreenArea = zeros(1, focusedN); % 
    pGreenLevel = zeros(1, focusedN); % 
    
    for sP = 1:focusedN
        % use idxBlue find which slice to use
        % green channels are in 2,4,6,8th... slice
        greenPickNum = idxGreen(sP)*dataChannels;
        pickedSliceGreen = dataImage{greenPickNum, 1}; % choose focused slice
        
        % level and segment
        levelGreen = adaptthresh(pickedSliceGreen, thresh_sensitivity); % use sensitivity=0.4 
        bw_green = imbinarize(pickedSliceGreen, levelGreen);
%         figure, imshow(bw);
        bw2_green = bwareaopen(bw_green, 3); % remove random noise (single pixel)
%         figure, imshow(bw2);

        % get corresponding blue channel, threshold it  
        pickedGreenBlue = dataImage{greenPickNum-1, 1}; % choose focused slice
        % level and segment
        level_greenBlue = adaptthresh(pickedGreenBlue, thresh_sensitivity); % use sensitivity=0.4   
        bw_greenBlue = imbinarize(pickedGreenBlue, level_greenBlue);
%         figure, imshow(bw);
        bw2_greenBlue = bwareaopen(bw_greenBlue, 3); % remove random noise (single pixel)
%         figure, imshow(bw2);
       % remove background, use uniform averaged level
       greenBlueBkg = mean(level_greenBlue(:))*(2^16);
       pick_greenBlueRvBkg = imsubtract(pickedGreenBlue, greenBlueBkg); % unit16
       
       % calculate results
        pBlueOnGreen(1, sP) = sum(pick_greenBlueRvBkg(bw2_green)); % count the Red channel using green channel mask
        pBlueTotal (1, sP) = sum(pick_greenBlueRvBkg(bw2_greenBlue)); 
        pGreenArea(1, sP) = sum(bw2_green(:));
        pGreenLevel(1, sP) = mean(levelGreen(:));        
    end    
    
    % write to final data table
    dataCountTable(iData, 13:16) = pBlueOnGreen;
    dataCountTable(iData, 19:22) = pBlueTotal;
    dataCountTable(iData, 25:28) = pGreenArea;
    dataCountTable(iData, 31:34) = pGreenLevel;
    
    dataCountTable(iData, 1:4) = pBlueOnGreen ./ pGreenArea; %   blue count on GFP mask
    dataCountTable(iData, 7:10) = pBlueOnGreen ./ pBlueTotal; % (blue count on GFP mask) out of total
    
    disp(iData);
    
end

%% write file table, adding annotations
topLineTxt = {'dataName', 'recruitInt', 'recruitInt', 'recruitInt', 'recruitInt', ...
              'n', 'n',   'recruitPrct',  'recruitPrct',  'recruitPrct',  'recruitPrct', ...
              'n', 'n',   'blueOnGrn',  'blueOnGrn',  'blueOnGrn',  'blueOnGrn', ...
              'n', 'n',   'blueTotal',  'blueTotal',  'blueTotal',  'blueTotal', ...
              'n', 'n',   'greenArea',  'greenArea',  'greenArea',  'greenArea', ...
              'n', 'n',   'greenLevel',  'greenLevel',  'greenLevel',  'greenLevel', ...
              'n', 'n',   'blueCount',  'blueCount',  'blueCount',  'blueCount', ...
              'n', 'n',   'blueArea',  'blueArea',  'blueArea',  'blueArea', ...
              'n', 'n',   'blueLevel',  'blueLevel',  'blueLevel',  'blueLevel', ...
              };
          
% dataLabel, first column in final table         
dataLabel = cell(dataN, 1);
 for i = 1:dataN
     dataLabel{i, 1} = dataList(i).name;
 end
 
% cat together for final table
finalTable = cat(2, dataLabel, num2cell(dataCountTable));
finalTable = cat(1, topLineTxt, finalTable);

% write to exel
xlswrite(resultFileName, finalTable);







