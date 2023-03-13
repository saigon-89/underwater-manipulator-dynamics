function manip_animation(q,q_vect,Tr,dt)
n = numel(q);
figure
plot3(0, 0, 0)
xlabel('X'), ylabel('Y'), zlabel('Z')
grid on
axis equal
xlim([-1 1])
ylim([-1 1])
zlim([-1 2])
hold on
for i = 1:length(q_vect)
    delete(get(gca, 'Children'));
    x = 0; y = 0; z = 0;
    for j = 1:n
        tmp = subs(Tr{j}(1:3,4), q, q_vect(i,:)');
        x = [x; tmp(1)]; y = [y; tmp(2)]; z = [z; tmp(3)];
    end
    plot3(x, y, z, 'r*'), hold on
    plot3(x, y, z, 'r-'), hold on
    pause(dt/10);    
end

end