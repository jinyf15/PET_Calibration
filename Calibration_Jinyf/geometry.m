size=20.9;
dx = size/32;
dy = size/32;
x=-size/2:dx:size/2;
y=-size/2:dy:size/2;
plot([x(1) x(1)],[y(1) y(end)],'k');
hold on;
plot([x(1) x(end)],[y(1) y(1)],'k');
plot([x(end) x(end)],[y(1) y(end)],'k');
plot([x(1) x(end)],[y(end) y(end)],'k');
for i = 2:length(x)-1
    plot([x(i) x(i)],[y(1) y(end)],'r');
    plot([x(1) x(end)],[y(i) y(i)],'r');
end
axis equal
axis([x(1) x(end) y(1) y(end)])
hold off