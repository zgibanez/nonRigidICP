function LB_test(LB)

nCol = size(LB,2);
LB = LB(1:3:end,:);
LB = 100*(abs(LB)-abs(repmat(LB(:,1),1,nCol)))./abs(repmat(LB(:,1),1,nCol));
ctr = [0:500:10000];
figure();
for i=1:nCol-1
    t = strcat('It ',num2str(i));
    subplot(1,nCol,i);
    histogram(LB(:,i+1),ctr);
    title(t);
    axis([0 10000 0 5000]);
end




end

