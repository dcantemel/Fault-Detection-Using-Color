function idx = simIdx(bw,bwT)

[rows, cols] = size(bw);
dist = zeros(rows, cols);
R = 3;
n = 0;
for i = 1+R:rows-R
    for j = 1+R:cols-R
        if bwT(i,j) ~=0
            for l = 0:R
                tmp = bw(i-l:i+l,j-l:j+l);
                if sum(tmp(:))>0
                    dist(i,j) = l-0.5;
                    n = n+1;
                    break;
                end
                dist(i,j) = R-0.5;
            end
        end
    end
end
idx = sum(dist(:))/n;

% function idx = simIdx(bw,bwT)
% 
% [rows, cols] = size(bw);
% dist = zeros(rows, cols);
% R = 3;
% n = 0;
% for i = 1+R:rows-R
%     for j = 1+R:cols-R
%         if bw(i,j) ~=0
%             for l = 0:R
%                 tmp = bwT(i-l:i+l,j-l:j+l);
%                 if sum(tmp(:))>0
%                     dist(i,j) = l-0.5;
%                     n = n+1;
%                     break;
%                 end
% %                 dist(i,j) = R-0.5;
%             end
%         end
%     end
% end
% idx = sum(dist(:))/n;