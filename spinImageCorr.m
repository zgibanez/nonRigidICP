function C = spinImageCorr(img1,img2,lambda)
%Input:
% -img1, img2: two square matrices of the same dimensions.


%First all the bins that have the value 0 for either images are discarded
%since it would inflate the correlation
idx = find(img1 & img2);

% %Normalize values
% img1(idx) = img1(idx)/max(img1(idx));
% img2(idx) = img2(idx)/max(img2(idx));

R = corr2(img1,img2);
C = atanh(R)^2-lambda*(1/(numel(idx)-3));


% N = size(img1,1);
% 
% corr = N*sum(img1(idx).*img2(idx))-sum(img1(idx))*sum(img2(idx)) / sqrt((N*sum(img1(idx).^2)-sum(img1(idx))^2)*(N*sum(img2(idx).^2)-sum(img2(idx))^2)); 

end