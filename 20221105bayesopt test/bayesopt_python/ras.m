function result=ras(xx,yy)
    result=-(20 + xx.^2 + yy.^2 - 10*(cos(2*pi*xx) + cos(2*pi*yy)));
end
