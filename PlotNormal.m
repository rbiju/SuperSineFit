function PlotNormal(filename)
load(filename); %#ok<LOAD>
p = transpose(periodValues);
rA = transpose(radAccelVals);
tA = transpose(tanAccelVals);
lA = transpose(longAmpVals);

pd1 = fitdist(p, 'Normal');
pd2 = fitdist(rA, 'Normal');
pd3 = fitdist(tA, 'Normal');
%pd4 is with radial resdiuals
pd5 = fitdist(lA, 'Normal');

%plot fits obtained from 100% sample superimposed on data
figure(1);
answ = fminunc(fh, x00, options);
radfit = (answ(4)*sin(2*pi.*(xdata/answ(1))+answ(3))+answ(5)) + (answ(6)*sin(2*pi.*xdata/answ(2)+answ(7)));
tanfit = (answ(4)*sin(2*pi.*xdata/answ(1)+(answ(3)-(pi/2)))+answ(8)) + (answ(9)*sin(2*pi.*xdata/answ(2)+answ(10)));

sz = 7;
red = [1, 0.6, 1];
scatter(xdata,radData,sz,red,'filled');

hold on
blue = [0.3010, 0.7450, 0.9330];
scatter(xdata,tanData,sz, blue,'filled');
plot(xdata,radfit,'k',xdata,tanfit,'b');

hold off

title('Least Square Fits of Radial and Tangential Accelerations on Six Flags Carousel');
xlabel('Time(s)');
ylabel('Acceleration(m/s^2)');
legend({'Radial Data', 'Tangential Data', 'Radial Fit','Tangential Fit'});
saveas(gcf, [pwd '/RotatedFitPDFs/SimFit.pdf']);

%plot histogram of period values with gaussian fit superimposed 
figure(3);
histfit(p,25);
ylabel({'Occurrence'});
xlabel({'Period (s)'});
title({'Histogram of Period Values'});
stats = sprintf('\x03bc: %f \n\x03c3: %f \nPeriod(1\x03c3 Confidence): \n(%.2f +/- %.2f) s',pd1.mu, pd1.sigma, pd1.mu, pd1.sigma);
annotation('textbox',[.6 .5 .3 .3],'String',stats,'EdgeColor','none');
ylim([0 100]);
saveas(gcf, [pwd '/RotatedFitPDFs/PeriodHist.pdf']);

%plot histogram of Radial Acceleration values with gaussian fit superimposed 
figure(4);
histfit(rA,25);
ylabel({'Occurrence'});
xlabel({'Radial Acceleration (m/s^2)'});
title({'Histogram of Radial Acceleration Values'});
stats = sprintf('\x03bc: %f \n\x03c3: %f \nRadial Acceleration(1\x03c3 Confidence): \n(%.4f +/- %.4f) m/s^2',pd2.mu, pd2.sigma, pd2.mu, pd2.sigma);
disp(pd2.mu);
disp(pd2.sigma);
annotation('textbox',[.6 .5 .3 .3],'String',stats,'EdgeColor','none');
ylim([0 100]);
saveas(gcf, [pwd '/RotatedFitPDFs/RadialAccelerationHist.pdf']);

%plot histogram of Tangential Acceleration values with gaussian fit superimposed 
figure(5);
histfit(tA,25);
ylabel({'Occurrence'});
xlabel({'Tangential Acceleration (m/s^2)'});
title({'Histogram of Tangential Acceleration Values'});
stats = sprintf('\x03bc: %f \n\x03c3: %f \nTangential Acceleration(1\x03c3 Confidence): \n(%.4f +/- %.4f) m/s^2',pd3.mu, pd3.sigma, pd3.mu, pd3.sigma);
disp(pd3.mu);
disp(pd3.sigma);
annotation('textbox',[.6 .5 .3 .3],'String',stats,'EdgeColor','none');
ylim([0 100]);
saveas(gcf, [pwd '/RotatedFitPDFs/TangentialAccelerationHist.pdf']);

%plot radial residuals
figure(6);
resRadVec = radData - radfit;
scatter(xdata, resRadVec, 'x');
ylabel({'Residuals (m/s^2)'});
xlabel({'Time (s)'});
title({'Radial Residuals'});
orient(gcf,'landscape');
print(gcf,[pwd '/RotatedFitPDFs/RadialResiduals.pdf'],'-dpdf', '-fillpage');

%plot tangential residuals
figure(7);
resTanVec = tanData - tanfit;
scatter(xdata, resTanVec, 'x');
ylabel({'Residuals (m/s^2)'});
xlabel({'Time (s)'});
title({'Tangential Residuals'});
orient(gcf,'landscape');
print(gcf,[pwd '/RotatedFitPDFs/TangentialResiduals.pdf'],'-dpdf', '-fillpage');

%plot total residuals
figure(8);
resTotVec = resTanVec + resRadVec;
scatter(xdata, resTotVec, 'x');
ylabel({'Residuals (m/s^2)'});
xlabel({'Time (s)'});
title({'Total Residuals'});
orient(gcf,'landscape');
print(gcf,[pwd '/RotatedFitPDFs/TotalResiduals.pdf'],'-dpdf', '-fillpage');

%plot radial residual histogram
figure(9);
pd4 = fitdist(resRadVec, 'Normal');
histfit(resRadVec, 25);
ylabel({'Occurrence'});
xlabel({'Residual (m/s^2)'});
title({'Radial Residual Histogram'});
stats = sprintf('\x03bc: %.9f \n\x03c3: %f', pd4.mu, pd4.sigma);
annotation('textbox',[.6 .5 .3 .3],'String',stats,'EdgeColor','none');
disp(pd4.mu);
disp(pd4.sigma);
saveas(gcf, [pwd '/RotatedFitPDFs/RadialResidualHist.pdf']);

%%plot normal distribution of longamp values obtained from ChiSquared()
%print mu, sigma, and scaled uncertainty
figure(10);
histfit(lA,25);
ylabel({'Occurrence'});
xlabel({'Amplitude (m/s^2)'});
title({'Histogram of Amplitude Values'});
stats = sprintf('\x03bc: %f \n\x03c3: %f \nAmplitude(1\x03c3 Confidence): \n(%.5f +/- %.5f)',pd5.mu, pd5.sigma, pd5.mu, pd5.sigma);
annotation('textbox',[.6 .5 .3 .3],'String',stats,'EdgeColor','none');
ylim([0 100]);
saveas(gcf, [pwd '/RotatedFitPDFs/LongAmpHist.pdf']);

%summary txt file
sum_file = [pwd '/RotatedFitPDFs/summary.txt'];
wh = fopen(sum_file, 'w');
fprintf(wh, 'Bootstrap of 100%% of points with locked large period, large amplitude, large period, small period on rotated data: \n');
fprintf(wh, 'Fitting to: \n(x(4) * sin((2*pi*t)/x(1) + x(3)) + x(5)) + (x(6) * sin((2*pi*t)/x(2) + x(7))) \n');
fprintf(wh, '(x(4) * sin((2*pi*t)/x(1) + (x(3)-pi/2)) + x(8)) + (x(9) * sin((2*pi*t)/x(2) + x(10))) \n');
fprintf(wh, 'Initial Guesses --> (x): \n');
guesses = sprintf('%f \n', x00);
fprintf(wh, '%s', guesses);
fprintf(wh, 'Period: \n(%.2f +/- %.2f) s \n',pd1.mu, pd1.sigma);
fprintf(wh, 'Radial Acceleration: \n(%.4f +/- %.4f) m/s^2 \n',pd2.mu, pd2.sigma);
fprintf(wh, 'Tangential Acceleration: \n(%.4f +/- %.4f) m/s^2 \n',pd3.mu, pd3.sigma);
fprintf(wh, 'Amplitude: \n(%.5f +/- %.5f) \n',pd5.mu, pd5.sigma);
fclose(wh);
end
