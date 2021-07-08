%% optic disc image and loc 8
opticDisc = imread("opticDisc_cropped.BMP");
loc8 = imcrop(imread("C:\Users\orica\Downloads\E0post_LOC8.BMP"));
%%
loc8_2 = loc8;
loc8_gray = rgb2gray(loc8);
loc8_gray_double = double(loc8_gray);loc8_gray_double(80:end,:) = loc8_gray_double(80:end,:) + 0.6*loc8_gray_double(1:39,:);loc8_gray_double = loc8_gray_double + fliplr(loc8_gray_double);
[row,col] = find(loc8_gray_double>mean(loc8_gray_double(:))+1.5*std(loc8_gray_double(:)));
loc8_binary = zeros(size(loc8_gray));
for i = 1:length(row)
    loc8_2(row(i),col(i),:) = [0,255,0];
end
tmp = uint8(zeros(size(opticDisc)));tmp(8:8+117,14:14+164,:) = loc8_2;
loc8_3 = 0.8*opticDisc+0.2*tmp;
imwrite(loc8_3,"loc8_stimuli.bmp");
%% 
[centersBright, radiiBright] = imfindcircles(loc8_2,[1 20],'ObjectPolarity','bright');
figure;imshow(loc8_2);
viscircles(centersBright, radiiBright,'Color','b');

%% loc 9 moving bars
retina = imcrop(imrotate(imread("C:\Users\orica\Downloads\retina_comparison_article_ori.BMP"),20));
%%
movingBars = imcrop(imrotate(imread("C:\Users\orica\Downloads\moving_bars_comparison_article_ori.BMP"),20));
%%
movingBars_2 = movingBars;
movingBars_gray = rgb2gray(movingBars);
movingBars_gray_double = double(movingBars_gray);
movingBars_gray_double = movingBars_gray_double + circshift(movingBars_gray_double,50,2) + circshift(movingBars_gray_double,-50,2);
[row,col] = find(movingBars_gray_double>mean(movingBars_gray_double(:))+0.8*std(movingBars_gray_double(:)));
for i = 1:length(row)
    movingBars_2(row(i),col(i),:) = [0,255,0];
end
tmp = uint8(zeros(size(retina2)));tmp(:,15:15+148,:) = movingBars_2;
movingBars_3 = 0.8*retina2+0.4*tmp;
imwrite(movingBars_3,"loc9_stimuli.bmp");
