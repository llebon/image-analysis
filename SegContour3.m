% Segmentation routine written by Michael Elowitz and Yaron Antebi
function imseg = SegContour3(im)
im = double(im);

% parameters:

if nargin < 2,
    hminscalar = 400;
end

% First, identify the cell regions:
initial_threshold=mode(im(im>0));
noiselevel=fminsearch(@(t) countblob(im>t), initial_threshold);
localthreshold=1.5*noiselevel;

imseg=im>localthreshold;
imsegc = imclose(imseg,strel('disk',3));

% Second, use contour dilations to grow the centers out to boundaries:

disk1=strel('disk',1);

% normalize the peaks to the same level
imex=imextendedmax(im,2*median(double(im(:))));
mask=(watershed(bwdist(imex)))>0;
L=imex.*im;
segarea=0;
for i=1:50
    L=(imdilate(L,disk1).*mask);
end
Lf=(imdilate(L,disk1));
Lf(Lf==0)=max(max(Lf));
imn=imtophat(im./imfilter(Lf,fspecial('gauss',5,5)),strel('disk',50));

%imex=imextendedmax(imfilter(im,fspecial('gauss',3,3)),2*median(double(im(:))));
%expand to fill contours
imex=imdilate(imex,strel('disk',2));
for th=linspace(1,0,10)
    imexarea=0;
    while imexarea<sum(sum(imex))
        imexarea=sum(sum(imex));
        imex=bwmorph(imex,'thicken').*(imn>th).*imsegc | imex;
    end
end

imseg=imex.*imsegc;
imseg=imfill(imdilate(bwlabel(imex),disk1));
imsegperim=(imseg-imdilate(imseg,disk1))<0;
imseg(imsegperim==1)=0;

function num = countblob (M)
    [~,num]=bwlabel(M);
    num=-num;
