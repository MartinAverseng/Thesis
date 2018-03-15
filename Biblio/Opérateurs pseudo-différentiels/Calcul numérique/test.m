%% Aspect sur un triangle

clear all
close all
tri = rand(2,3); %trois pts au hasard
tri = [tri tri(:,1)];
sides = tri(:,[2 3 4])-tri(:,[1 2 3]);

lengths = sqrt(sum(sides.^2,1));
Perim = sum(lengths);
plot(tri(1,:) ,tri(2,:));

N = 3000;
t0 = fix(N*rand(1));
t = linspace(0,1,N);

for i = 1:length(sides)
    if i == 1
        thisSideIndexes = t*Perim<=lengths(i);
    else
        thisSideIndexes = and(sum(lengths(1:i-1))<t*Perim,Perim*t<=sum(lengths(1:i)));
    end
      
       tSide = t(thisSideIndexes);
       y1(thisSideIndexes)= tri(1,i) + (tri(1,i+1)-tri(1,i))*(tSide-sum(lengths(1:i-1))/Perim)/lengths(i)*Perim;
       y2(thisSideIndexes)= tri(2,i) + (tri(2,i+1)-tri(2,i))*(tSide-sum(lengths(1:i-1))/Perim)/lengths(i)*Perim;
end

hold on
plot(y1,y2,'r--');
x = [y1(t0),y2(t0)];
dist = sqrt((x(1)-y1).^2+(x(2)-y2).^2);
plot(x(1),x(2),'*');
figure;
plot(t,log(dist));
hold on;


for k = 1:100;
    an(k) = mean((dist).*exp(2*1i*k*pi*t));   
end
    
figure
plot(log(abs(an)));
    
    