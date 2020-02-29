function error = funfit4(guesses, xdata, ydata)
y1 = guesses(1) + guesses(2).*sin(((2.*pi.*xdata)./guesses(3)) + guesses(4)) + guesses(5).*sin(((2.*pi.*xdata)./guesses(6)) + guesses(7));
%    Short, Long
error = sum((ydata-y1).^2);
end