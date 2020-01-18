%{ 
ELEC 4700: Assignment 1
William Fraser
101001393
%}

numParticles = input("Please input the number of particles: ");
kb = 1.38064852e-23;
m0 = 9.11e-31;
m = 0.26*m0;
T = 300;
kqq = (1e18)*(1.609e-19)^2;

vth = sqrt((kb*T)/m);

grid1 = zeros(1,10000);
initParticles = zeros(numParticles,3);
initParticles(:,2) = vth*1e-5;

for i = 1:length(initParticles(:,1))
    initParticles(i,1) = round(rand*1000);
    randomDirection = (-1)^(round(rand));
    initParticles(i,2) = initParticles(i,2)*randomDirection;
end

uniquePositions = unique(initParticles(:,1));
keepPoints = zeros(length(initParticles(:,1)),1);
for i = 1:length(uniquePositions)
    temp = initParticles(:,1) == uniquePositions(i);
    found = 0;
    currentIndex = 1;
    while(found==0)
        if(temp(currentIndex)==1)
            found = 1;
        else
            currentIndex = currentIndex + 1;
        end
    end
    keepPoints(currentIndex) = 1;
end

particles(:,1) = nonzeros(initParticles(:,1).*keepPoints);
particles(:,2) = nonzeros(initParticles(:,2).*keepPoints);
particles(:,3) = 0;

particles = sortrows(particles,1);

for i = 1:length(particles(:,1))
    position = particles(i,1);
    grid1(position+1) = 1;
end

spy(grid1,10)
force = zeros(length(particles(:,1)),1);

for time = 1:1000
    tempParticles = particles;
    grid2 = zeros(1,1000);
    posVector = particles(:,1);
    left = circshift(posVector,1);
    right = circshift(posVector,-1);
    stopPoint = length(particles(:,1)) - 1;
    force(2:stopPoint) = (kqq./((abs(particles(2:stopPoint,1)-left(2:stopPoint))).^2)) - (kqq./((abs(particles(2:stopPoint,1)-right(2:stopPoint))).^2));
    force(1) = -kqq/(abs(particles(1,1)-right(1))^2);
    force(length(particles(:,1))) = kqq/(abs(particles(length(particles(:,1)),1)-left(length(particles(:,1))))^2);
    force(isnan(force)) = 0;
    particles(:,3) = (force./m0)*1e-10;
    particles(:,2) = particles(:,2) + particles(:,3);
    tempParticles(:,1) = particles(:,1) + particles(:,2);
    goodPoints1 = tempParticles(:,1) > 1;
    badPoints1 = tempParticles(:,1) < 1;
    goodPoints2 = tempParticles(:,1) < 1000;
    badPoints2 = tempParticles(:,1) > 1000;
    particles(:,1) = (tempParticles(:,1).*goodPoints1.*goodPoints2) + (badPoints1*1) + (badPoints2.*1000);
    particles = sortrows(particles,1);
    for i = 1:length(particles(:,1))
        position = round(particles(i,1));
        if(position==0)
            position = position + 1;
        end
        if(~isnan(position))
            grid2(position) = 1;
        end
    end
    grid1 = grid2;
    clear grid2;
    clear tempParticles
    spy(grid1,10)
    pause(0.1)
end
    
    
