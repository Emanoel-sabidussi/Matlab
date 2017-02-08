close all
clear all

tic
im_or = imread('100.tif');
im=255*im2double(im_or);

gradIm = gradient(double(im));
figure
imshow(gradIm)
indx = find(abs(gradIm)>mean(mean(abs(gradIm)))*2.2263);

gradIm2 = zeros(size(im_or,1),size(im_or,2));

reshGrad = reshape(gradIm2, size(im_or,1)*size(im_or,2), 1);
reshGrad(indx) = 1;

imGrad = reshape(reshGrad,size(im_or,1),size(im_or,2));
figure
imshow(imGrad)

se = strel('disk',3);
dilateBW = imdilate(imGrad, se);
se = strel('disk',4);
erodedBW = imerode(dilateBW, se);
figure
imshow(erodedBW)

I0 = imfill(erodedBW,'holes');
I0 = imerode(I0, se);
I0 = imdilate(I0, se);



I=double(I0);

s=1;
I = gaussianBlur(I,s);

%parameters
alpha=0.2;
mu=0.2;
iter=80;

% run GVF
tic;
[u,v] = GVF(I, alpha, mu, iter);
t=toc;

fprintf('Computing GVF uses %f seconds \n',t);

% visualize result
imshow(im_or);
hold on;
quiver(u,v);
axis ij off;

%im = imadjust(im,[0.1 0.7],[]);

toc
%%
ticindex_beg = [1:size(im,2)];
index_end = [1:size(im,2)];

im2 = im;

med = 10;
coeff24hMA = ones(1, med)/med;

line_ts = double(im);
md_int_line = mean(line_ts);
avg_24_line = filter(coeff24hMA, 1, line_ts);

for i=2:size(im,2)
   
    idx = find(avg_24_line(:,i)>1*md_int_line(:,i));
    
    if isempty(idx)
        index_beg(i) = index_beg(i-1);
    else
        index_beg(i) = idx(1);
    end
    
end
toc

%%

im2 = wiener2(im,[25 25]);

mid_point = index_beg(round (size(im,2)/2)-2);
lineL = 1:1:size(im,2);

ang = pi/2;
for i=1:90
    
   XL = (mid_point+(ang*(size(im,2)/2)))-(0.7848*size(im,2)) + lineL*cos(ang);
   ang = ang-0.001;
   
   for j =1:(size(im,2))
        image_lin1(j) = im(round(abs(XL(j))+1),lineL(j));
   end
   linefin1(:,i) = (round(abs(XL)));
   mean1(i) = mean(image_lin1);
end



ang = pi/2;
for i=90:180
   
   XL = (mid_point+(ang*(size(im,2)/2)))-(0.7848*size(im,2)) + lineL*cos(ang);    
   ang = ang+0.001;
   
   for j =1:(size(im,2))
        image_lin2(j) = im(round(abs(XL(j))+1),lineL(j));
   end
   linefin1(:,i) = XL;
   mean1(i) = mean(image_lin2);
   
end

line_indx1 = max(mean1);

fin_indx = find(mean1==(max(mean1)));   
 
im3 = im;

   for j =1:(size(im,2))
        im3(    round(linefin1(j,(fin_indx)))-4 : round(linefin1(j,(fin_indx)))+4      ,lineL(j)) = im(    round(linefin1(j,(fin_indx)))-4 : round(linefin1(j,(fin_indx)))+4      ,lineL(j)) - im(    round(linefin1(j,(fin_indx)))-4 : round(linefin1(j,(fin_indx)))+4      ,lineL(j));
   end

%hold on
%plot(linefin1(:,fin_indx))
%

%%

tic
line_ts1 = double(I0);
md_int_line1 = mean(line_ts1);
avg_24_line1 = filter(coeff24hMA, 1, line_ts1);

for i=2:size(im,2)
   
    idx1 = find(avg_24_line1(:,i)>1.1*md_int_line1(:,i)); 
    idx2 = find(avg_24_line1(:,i)>1.7*md_int_line1(:,i));
    if isempty(idx1)
        index_beg(i) = index_beg(i-1);
    else
        index_beg(i) = idx1(1);
    end
    
    if isempty(idx2)
        index_end(i) = index_end(i-1);
    else
        index_end(i) = idx2(end);
    end
    
end
toc
%%
tic
FILT = load('FILT.mat');
FILT1 = load('FILT1.mat');
FILT = FILT.FILT;
FILT1 = FILT1.FILT1;


index_beg_filt = filtfilt(FILT,index_beg);
index_end_filt = filtfilt(FILT1,index_end);
toc

%%
tic
[index_end_filt_env_hi, index_end_filt_env_lo] = envelope(index_end_filt,25,'peak');


figure
imshow(im_or)
hold on

i = linspace(1,size(im,2),size(im,2));

plot(i,index_beg_filt,'r.','MarkerSize',5);
plot(i,index_end_filt_env_lo,'b.','MarkerSize',5);

toc


%%

length = abs(index_beg_filt(1:1000) - index_end_filt_env_lo(1:1000));
length = length*720/size(im,2);
length = mean(length)
