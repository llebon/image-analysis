function [segmented] = segImagesDAPI(s_im ,r_im, y_im)

% This function takes in three .tif images from three fluorescence channels.
% s_im is the segmentation image that is used to identify cell bodies.  
% In this experiment, we used CellMask Blue cytoplasmic stain to label cell bodies, and used the DAPI channel to image this stain.
% The other two images are the fluorescence images (YFP and RFP channels) from which cell fluroscence data is extracted.  
% The properties for each cell in the image are stored in the "segmented" data structure.

% Read in the images

seg_im = imread(s_im); % DAPI channel is the segemntation color in this experiment
rfp_im = imread(r_im);
yfp_im = imread(y_im);

[m, n] = size(seg_im); % Calculate image size

% Initialize the output data structure

segmented = struct('coords', {}, 'area', {}, 'rect', {}, 'bw', {},...
    'rfp_median', {}, 'yfp_median', {}, 'seg_median', {},  ...
    'rfp_mean', {}, 'yfp_mean', {}, 'seg_mean', {},  ...
    'rfp_tot', {}, 'yfp_tot', {}, 'seg_tot', {},  ...
    'yfp_peak', {}, ...
    'rfp_bg', {}, 'yfp_bg', {}, 'seg_bg', {}) ;

% Segment the DAPI image.

LN = SegContour3(seg_im); % Michael and Yaron's segmentation
LN = ourclearborder(LN); % Omit the cells on the edges of the image

% Collect the statistics for each cell.

Stats = regionprops(LN, 'Area', 'Centroid', 'BoundingBox');

% For each cell in the image, collect the area, centroid, bounding box, and fluorescence data

for i = 1:max(max(LN));
    imCellD = double(seg_im(LN==i));
    imCellR = double(rfp_im(LN==i));
    imCellY = double(yfp_im(LN==i));
    
    % Evaluate the centroid position of the cell
    
    centroid = Stats(i).Centroid;
    
    % Evaluate the peak YFP value in the cell (Top 1% of pixel values in
    % the YFP channel)
    
    sortYFP = sort(imCellY, 'descend');
    peakYFP = mean(sortYFP(1:(round(.01*Stats(i).Area))));
    
    % Evaluate the bounding box of the cell
    
    rect = Stats(i).BoundingBox;
    Lcrop = imcrop(LN, rect);
    segmented(i).bw = Lcrop;
    
    % Save the area, centroid, and bounding box of the cell to the output
    % structure
    
    segmented(i).area = Stats(i).Area;
    segmented(i).coords = Stats(i).Centroid;
    segmented(i).rect = Stats(i).BoundingBox;
    
    % Calculate local fluorescence backgrounds - median of bottom 25% of pixels in the
    % bounding box
    
    rectD = imcrop(seg_im, rect);
    rectR = imcrop(rfp_im, rect);
    rectY = imcrop(yfp_im, rect);
    
    rectD = double(rectD);
    rectR = double(rectR);
    rectY = double(rectY);
    
    sortedD = sort(rectD);
    sortedR = sort(rectR);
    sortedY = sort(rectY);
    
    minImD = median(sortedD(1:round(.25*length(sortedD))));
    minImR = median(sortedR(1:round(.25*length(sortedR))));
    minImY = median(sortedY(1:round(.25*length(sortedY))));
    
    
    % Evaluate the background values for each channel
    
    segmented(i).seg_bg = minImD;
    segmented(i).rfp_bg = minImR;
    segmented(i).yfp_bg = minImY;
    
    % Evaluate the median fluorescence for each channel
    
    fim_medD = median(imCellD(:) -  minImD);
    fim_medR = median(imCellR(:) - minImR);
    fim_medY = median(imCellY(:) - minImY);
    
    % Evaluate the mean fluorescence for each channel
    
    fim_meanD = mean(imCellD(:) - minImD);
    fim_meanR = mean(imCellR(:) - minImR);
    fim_meanY = mean(imCellY(:) - minImY);
    
    % Evaluate the total fluorescence for each channel
    
    seg_total = sum(imCellD(:)) - minImD*length(imCellD(:));
    rfp_total = sum(imCellR(:)) - minImR*length(imCellR(:));
    yfp_total = sum(imCellY(:)) - minImY*length(imCellY(:));
    
    % Evaluate the peak YFP value
    
    segmented(i).yfp_peak = peakYFP;
    
    % Save all the fluorescence values for this cell to the output
    % structure
    
    segmented(i).seg_median = fim_medD;
    segmented(i).rfp_median = fim_medR;
    segmented(i).yfp_median = fim_medY;
    
    segmented(i).seg_mean = fim_meanD;
    segmented(i).rfp_mean = fim_meanR;
    segmented(i).yfp_mean = fim_meanY;
    
    segmented(i).seg_tot = seg_total;
    segmented(i).rfp_tot = rfp_total;
    segmented(i).yfp_tot = yfp_total;
    

    clear fim_medR fim_medY fim_meanR fim_meanY rfp_total yfp_total seg_total;
    clear fim_meanD fim_medD;
    clear rectD sortedD minImD;
    clear rectR rectY sortedR sortedY sortYFP peakYFP minImR minImY rect
    clear imCellR imCellY centroid peakYFP sortYFP imCellD;
end
