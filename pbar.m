function s = pbar(n,N)
%PBAR Creates a progressbar
stars = round(n/N*20);
dots = 20 - stars;
s = strcat('[',repmat('*',1,stars),repmat('·',1,dots),']');
end

