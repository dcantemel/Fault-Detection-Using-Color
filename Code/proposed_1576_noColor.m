% load 3D data sets
addpath functions
addpath Skeleton
load dataTime.mat   % Original seismic data
load sembTime.mat   % Corresponding semblance attribute

TimeRes = 113;      % Time Resolution(z axis)
NTrace = 151;       % Number of traces (x axis)
NInline = 300;      % Number of inline (y axis)

% The dimension of the test image
rows = NTrace;         
cols = NInline;

% Parameters in adapthisteq
NumTiles = [2, 2];
ClipLimit = 0.01;   
NBins = 100;        
Range = 'full';
Distribution='Rayleigh';
Alpha=0.4;

% Parameters used to fill the margins of images
rm = 3;
fillVal = 0;

% Parameters used to smooth the semblance map
filtSze = [2,2];
sigma = 10;
h = fspecial('gaussian', filtSze, sigma);

% Target section
Section = 44;

% Combine the neighboring sections and generate a RGB image
tmpSemb = semb(:,:,Section);
testSemb = tmpSemb/(max(tmpSemb(:)));

% figure,imshow(testSemb,'border', 'tight');

% Smooth the semblance map
testSemb1 = imfilter(testSemb, h);
% figure,imshow(testSembV1,'border','tight');

% Apply adaptive histogram equalization
testSemb2 = adapthisteq(testSemb1,  'NumTiles', NumTiles, 'ClipLimit', ClipLimit, 'NBins', NBins, 'Range', Range);
% figure,imshow(testSembV2,'border','tight');

idxV = testSemb2<0.5;
idxV = marginFilling(idxV, rm, fillVal);
figure,imshow(idxV, 'border','tight');

faultReg2 = idxV;

Data = data(:,:,Section);
a = 20000;
Data1 = Data;
tmp = find(abs(Data(:))>a);
Data1(tmp) = sign(Data(tmp)).*a;
Data2 = Data1./a;                           % Normalization from [-a, a] to [-1, 1]
data1 = Data2;
data2 = hilbert(data1);

semblance = semb(:,:,Section);
semblance = marginFilling(semblance, rm, max(semblance(:)));

semblance = semblance./(max(semblance(:)));
% figure, pcolor(semblance(4:end-4, 4:end-4));
% shading interp
% axis ij
% colorbar;
% colormap;
% axis equal tight
% set(gca, 'FontSize', 14);
% title('Discontinuity Map','fontsize', 18);
%%
r = 1;
% Simple skeletonization
bw0 = bwmorph(skeleton(faultReg2)>10,'skel',Inf);

% Weighted Average of discontinuity value
enve = abs(data2).^2;

vol = zeros(rows-2*r, cols-2*r, (2*r+1)^2);
weight = zeros(rows-2*r, cols-2*r, (2*r+1)^2);

for i = 1:2*r+1
    for j = 1:2*r+1
        vol(:,:,(2*r+1)*(i-1)+j) = semblance(i:rows-2*r+i-1,j:cols-2*r+j-1).*enve(i:rows-2*r+i-1,j:cols-2*r+j-1);
    end
end

for i = 1:2*r+1
    for j = 1:2*r+1
        weight(:,:,(2*r+1)*(i-1)+j) = enve(i:rows-2*r+i-1,j:cols-2*r+j-1);
    end
end

ave_discon = zeros(rows, cols);
ave_discon(r+1:rows-r, r+1:cols-r) = 1./sum(weight,3).*sum(vol,3);
bw1 = bwmorph(skeleton(faultReg2).*ave_discon>10,'skel',Inf);
figure,imshow(bw1,'border','tight');

% Eliminate the isolated points
bw2 = zeros(rows, cols);
for i = 1+r:rows-r
    for j = 1+r:cols-r
        if (sum(sum(bw1(i-1:i+1,j-1:j+1))) == 1) % indicate an isolated point
            bw2(i,j)=0;
        else
            bw2(i,j) = bw1(i,j);
        end
    end
end
% figure,imshow(bw2, [], 'border', 'tight');

% Eliminate the short line segments in bw2
bw3 = bw2;    % Final results after removing short line segments
flag = zeros(rows, cols);   % Flag matrix to indicate the points has been calculated with 1

r = 1;
pts = zeros(1,2);
lenThres = 8;
for i = 1+r:rows-r
    for j = 1+r:cols-r
        seed = bw2(i,j);
        idr = i;    % row index of seed
        idc = j;    % col index of seed
%         if i==62 && j == 54
%             ttt = 0;
%         end
        while(seed && (flag(idr, idc) == 0))
%             bw3(idr, idc) = 1;
            pts = [pts;idr idc];
            flag(idr, idc) = 1;
            
            % neighboring of seed
            bin = bw2(idr-r:idr+r, idc-r:idc+r);
            % show points in the neighboring of seed has been calculated or not
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
                seed = bw2(idr,idc);
            else
                seed = 0;
                len = length(pts(2:end, :)); 
                %_________________________________________
                bin1 = bw2(pts(2,1)-r:pts(2,1)+r, pts(2,2)-r:pts(2,2)+r);
                bin2 = bw2(pts(end,1)-r:pts(end,1)+r, pts(end,2)-r:pts(end,2)+r);
               
                [idx_bin1_r, idx_bin1_c] = find(bin1 == 1);
                idr_bin1 = pts(2,1)+idx_bin1_r-(r+1);
                idc_bin1 = pts(2,2)+idx_bin1_c-(r+1);
                idx_bin1_pts = (idc_bin1-1)*rows+idr_bin1;
                
                [idx_bin2_r, idx_bin2_c] = find(bin2 == 1);
                idr_bin2 = pts(end,1)+idx_bin2_r-(r+1);
                idc_bin2 = pts(end,2)+idx_bin2_c-(r+1);
                idx_bin2_pts = (idc_bin2-1)*rows+idr_bin2;
                 %_________________________________________
                 
                pts_idx = pts(2:end,1)+(pts(2:end,2)-1)*rows;
                
                 %_________________________________________
                 com1 = intersect(idx_bin1_pts, pts_idx);
                 len1 = length(idx_bin1_pts);
                 flag1 = (length(com1) == len1);
                
                 com2 = intersect(idx_bin2_pts, pts_idx);
                 len2 = length(idx_bin2_pts);
                 flag2 = (length(com2) == len2);
                 %_________________________________________
                if (len < lenThres) && ((flag1 == 1) || (flag2 == 1))
%                 if (len<=lenThres)
                    bw3(pts_idx) = 0;
                end
                pts = zeros(1,2);
            end
        end
    end
end


% figure,imshow(bw3,'border','tight');
% semb = Discon(:,:,Section);
% semb_idx = find(semb == 0);
% semb(semb_idx) = max(semb(:));

fimg =imfuse(bw3,semblance,'blend', 'scaling', 'joint');
figure, imshow(fimg,'border', 'tight');




