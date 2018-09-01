%% function to analyze the PDFs from the floats

% using ksdensity to estimate PDF
figure

for id =5:30
    
    maxsep = log10(max(dispersion.disp(id,:).^0.5)+1000);
    vecx = logspace(log10(3000),maxsep,20);
    
    [pdf, rpdf] = ksdensity(dispersion.disp(id,:).^0.5,vecx,'support','positive');
    
    %
    %bar(rpdf/1000,pdf,'edgecolor','k','facecolor','k')
    loglog(rpdf/1000,pdf)
    hold all
end
axis([1 1000 10^-10 1])
%         plot(r, dist_lung(:,timeid),'color','b')figure, clf,


%% raw PDFs 
figure

for id =10:20

[data_pdf,b] = hist(dispersion.disp(id,:).^0.5,vecx);

loglog(b/1000,data_pdf)
hold all
end