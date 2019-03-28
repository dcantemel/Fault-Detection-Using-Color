 % This is the first version for the implementation of Bo Zhang's paper
% (2013) "Semi-automated fault interpretation based on seismic attributes"
% In this paper, authors detected the faults in the time slices by
% borrowing the ideas of finger-vein pattern detection from medical image
% processing.
addpath functions
load sembTime.mat % load the semblance cube 
% addpath('F:\Z Wang\Research\Matlab\SegyMat\SegyMAT');
% load testTSlice402.mat
% load Semb_sqrt_old.mat
% [Data, SegyTraceHeaders,SegyHeader]=ReadSegy('F3.sgy');
% Description F3
% Inline: range: 200~650, step=2, NInline=226
% Xline: range: 700~1200, step=2,NTrace=251
% TimeRes: range: 400~1100, step=4, TimeRes=176

%%
% % Source code generates the discontinuity cube aligned on time slices
% % Read one seismic section from 3D Dataset
% TimeRes = 176;
% NTrace = 251;
% NInline = 226;
% % Data_sSlice: slices aligned on inline direction with the dimension of [TimeRes, NTrace, NInline]
% Data_sSlice = Seis2D23D(Data, TimeRes, NTrace, NInline);
% % Data_cSlice: slices aligned on crossline direction with the dimension of
% % [TimeRes, NInline, NTrace];
% Data_cSlice = Seis2D23DC(Data, TimeRes, NTrace, NInline);
% Semb_sSlice = semb3D(Data_sSlice);
% Semb_cSlice = semb3D(Data_cSlice);
% Semb_stSlice = Seis2Time(Discon_sSlice, TimeRes, NTrace, NInline);
% semb_ctSlice = Crs2Time(Discon_cSlice, TimeRes, NTrace, NInline);
% Semb = sqrt(Semb_stSlice.*Semb_ctSlice);
%%
% ----------Image-------------
% for i = 30:50
% figure, pcolor(Discon(:, :, i));
% shading interp
% axis ij
% colorbar
% end

%%
% Read one time section from 3D Dataset
% [rows, cols] = size(img);
% Data = img;
% 
% % Project seismic data region to [0,1]
% % Parameter 1: a
% % Firstly project seismic data region to [-a, a], a>0, here a = 20000
% % Then from [-a, a] to [-1, 1]
% a = 20000;
% Data1 = Data;
% tmp = find(abs(Data(:))>a);
% Data1(tmp) = sign(Data(tmp)).*a;
% Data2 = Data1./a;                           % Normalization from [-a, a] to [-1, 1]
% data = Data2;
% data1 = hilbert(data);

Section = 53;
% tmp = abs(log(Discon(:,:,Section)));
% 
% ind_inf = find(isinf(tmp(:)));
% tmp(ind_inf) = 0;

img = semb(:,:,Section);
% img = img./(max(img(:)));
[rows, cols] = size(img);

fvr = ones(rows, cols);
mask = find(img==0);

fvr(mask) = 0;

%%
% % Source code for finger-vein pattern detection
% img = im2double(imread('finger.png')); % Read the image
% img = imresize(img,0.5);               % Downscale image

% Get the valid region, this is a binary mask which indicates the region of 
% the finger. For quick testing it is possible to use something like:
% fvr = ones(size(img));
% The lee_region() function can be found here:
% http://www.mathworks.com/matlabcentral/fileexchange/35752-finger-region-localisation

% It's obvious that if we use smaller masks, we will get better results
% fvr = lee_region(img,4,40);    % Get finger region
% figure,imshow(fvr);

%% Extract veins using maximum curvature method
sigma = 6; % Parameter
v_max_curvature = miura_max_curvature(img,fvr,sigma);

% Binarise the vein image
md = median(v_max_curvature(v_max_curvature>0));
% v_max_curvature_bin = v_max_curvature > md; 
v = v_max_curvature > md;
% figure,imshow(v);

% img(mask)=max(img(:));
% figure,imshow(img);

v1 = bwmorph(v, 'skel', Inf);
% figure,imshow(v1,'border','tight');
r = 1;

v2 = zeros(rows, cols);
for i = 1+r:rows-r
    for j = 1+r:cols-r
        if (sum(sum(v1(i-1:i+1,j-1:j+1))) == 1) % indicate an isolated point
            v2(i,j)=0;
        else
            v2(i,j) = v1(i,j);
        end
    end
end
figure,imshow(v2, [], 'border', 'tight');

%%
% Eliminate the short line segments in bw2
v3 = v2;    % Final results after removing short line segments
flag = zeros(rows, cols);   % Flag matrix to indicate the points has been calculated with 1

r = 1;
pts = zeros(1,2);
lenThres = 8;
for i = 1+r:rows-r
    for j = 1+r:cols-r
        seed = v2(i,j);
        idr = i;    % row index of seed
        idc = j;    % col index of seed
        while(seed && (flag(idr, idc) == 0))
%             bw3(idr, idc) = 1;
            pts = [pts;idr idc];
            flag(idr, idc) = 1;
            
            % neighboring of seed
            % bin = [bw2(idr-r:idr-1,idc-r:idc+r), bw2(idr, idc+r:-1:idc+1), bw2(idr, idc-1:-1:idc-r), bw2(idr+1:idr+r, idc-r:idc+r)];
            bin = v2(idr-r:idr+r, idc-r:idc+r);
            % show points in the neighboring of seed has been calculated or not
            % sig = [flag(idr-r:idr-1,idc-r:idc+r), flag(idr, idc+r:-1:idc+1), flag(idr, idc-1:-1:idc-r), flag(idr+1:idr+r, idc-r:idc+r)];
            sig = flag(idr-r:idr+r, idc-r:idc+r);
            idx_v = find(bin == 1);
            idx_v(idx_v == 5) =[];
            idx_s = find(sig == 0);
            idx_s(idx_s == 5) =[];
            [com, IV, IS] = intersect(idx_v, idx_s);
            if ~isempty(com)
                idx = min(IV);
                idx_c = ceil(idx_v(idx)/(2*r+1));
                idx_r = idx_v(idx)-(idx_c-1)*(2*r+1);
                
                idr = idr+idx_r-(r+1);
                idc = idc+idx_c-(r+1);
                seed = v2(idr,idc);
            else
                seed = 0;
                len = length(pts(2:end, :)); 
                %_________________________________________
%                 bin1 = bw2(pts(2,1)-r:pts(2,1)+r, pts(2,2)-r:pts(2,2)+r);
%                 bin2 = bw2(pts(end,1)-r:pts(end,1)+r, pts(end,2)-r:pts(end,2)+r);
%                
%                 [idx_bin1_r, idx_bin1_c] = find(bin1 == 1);
%                 idr_bin1 = pts(2,1)+idx_bin1_r-(r+1);
%                 idc_bin1 = pts(2,2)+idx_bin1_c-(r+1);
%                 idx_bin1_pts = (idc_bin1-1)*rows+idr_bin1;
%                 
%                 [idx_bin2_r, idx_bin2_c] = find(bin2 == 1);
%                 idr_bin2 = pts(end,1)+idx_bin2_r-(r+1);
%                 idc_bin2 = pts(end,2)+idx_bin2_c-(r+1);
%                 idx_bin2_pts = (idc_bin2-1)*rows+idr_bin2;
                 %_________________________________________
                 
                pts_idx = pts(2:end,1)+(pts(2:end,2)-1)*rows;
                
                 %_________________________________________
%                  com1 = intersect(idx_bin1_pts, pts_idx);
%                  len1 = length(idx_bin1_pts);
%                  flag1 = (length(com1) == len1);
%                 
%                  com2 = intersect(idx_bin2_pts, pts_idx);
%                  len2 = length(idx_bin2_pts);
%                  flag2 = (length(com2) == len2);
                 %_________________________________________
%                 if (len < lenThres) && (flag1 == 1) && (flag2 == 1)
                if (len<lenThres)
                    v3(pts_idx) = 0;
                end
                pts = zeros(1,2);
            end
        end
    end
end

img(mask)=max(img(:));



%%
sel = strel('disk',1);
bw1 = imdilate(v3, sel);
bw2 = imerode(bw1, sel);
bw3 = bwmorph(bw2, 'skel', Inf);

fimg =imfuse(bw3,img,'blend', 'scaling', 'joint');
figure,imshow(fimg,'border', 'tight');




















