function [periodValues, radAccelVals, tanAccelVals, longAmpVals] = ChiSquared(filename)
%%function to generate bootstrap distribution for parameters of
%%simultaneous fit of two sine waves
tic %timer
periodValues = zeros(1,750); %arrays to store parameter of interest
radAccelVals = zeros(1,750);
tanAccelVals = zeros(1,750);
longAmpVals = zeros(1,750);

load(filename); %#ok<LOAD>
[~, c] = size(randNdxs); %#ok<*NODEF>
randRad = radData(randNdxs);
randTan = tanData(randNdxs);
global xdata; %#ok<NUSED>
for i = (1:c) %iterates through each pre-generated collection of random indexes(duplicates allowed)
    times = randNdxs(:,i); %#ok<IDISVAR>
    fun = @(x00) funFit(x00, randRad(:,i), randTan(:,i), times); %function handle to pass to fminunc()
    params = fminunc(fun, x00, options); %fminunc allowed to iterate 2000 times-- ideal amount unknown
    periodValues(i) = params(1);
    radAccelVals(i) = params(5);
    tanAccelVals(i) = params(8);
    longAmpVals(i) = params(4);
end
toc
end

function error = funFit(x, ydata1, ydata2, times)
perbig = x(1);
persmall = x(2);
phaseRad = x(3);
phaseTan = phaseRad - (pi/2); %locks both waves to be 90 degrees phase shifted
amp = x(4);
x1 = x(5:7);
x2 = x(8:10);
y1 = Voigt(x1, perbig, persmall, phaseRad, times,amp);
y2 = Voigt(x2, perbig, persmall, phaseTan, times,amp);
error = sum((ydata1-y1).^2 + (ydata2-y2).^2); %least squares error to minimize
end

function y = Voigt(x, perbig, persmall, phase, times,amp)
global xdata;
y = (amp * sin(2*pi.*xdata(times)./perbig + phase) + x(1)) + (x(2) * sin(2*pi.*xdata(times)./persmall+x(3))); %sum of sines function I am fitting to
end