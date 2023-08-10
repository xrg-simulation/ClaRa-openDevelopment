close all
clear all

%% assumptions

% no distance between collcetors 1-3 and 4-6
% no distance between evaporator and superheater field
% when the distortion reaches the beginning or the end of one irradiation discetization, the whole
% discretization is being shaded or insolated at once
% clouds have infinitive width

%% input data

% General
startingTimeDistortion=20000;
endTimeSim=40000;

%Distortion
orientation=180; %orientation of distortion in degrees north=0�, eastwards positive
lengthDistortion=300;
velocityDistortion=1;
deltaTime=20;
Irr=850;
deltaIrr=0.8*850;
nDistortion=1; %Number of Distortions in row
diffLengthDistortion=200; %Length of gap between distortions

%Evaporator
lengthTube=720;
discIrr=36;     %zones of different irradiation for evaporator
nString=6; %number of evaporator loops
widthEvapField=264; %width of solar evaporator field

%Superheater
nStringSH=3; %number of superheater loops
widthSHField=77; %width of solar superheater field
lengthTubeSH=240; %length of superheater tube
discIrrSH=8;     %zones of different irradiation for superheater


%% general calculations

deltaTimeDistortion=lengthDistortion/velocityDistortion;
velocityDistortionHor=abs(sind(orientation)*velocityDistortion);
velocityDistortionVert=abs(cosd(orientation)*velocityDistortion);
lengthTubeHor=abs(round(lengthTube/discIrr*sind(orientation)));
lengthTubeHorSH=abs(round(lengthTubeSH/discIrrSH*sind(orientation)));
lengthTubeHorSHPowerBlock=abs(round(lengthTubeSH/2*sind(orientation)));
startingTimeDistortionDummy2=startingTimeDistortion;

%% create table data evaporator

if (orientation <= 90) || (orientation >= 270)
    startingTimeDistortionDummy=startingTimeDistortion+round(widthEvapField*sind(90-orientation)/velocityDistortion);
else
    startingTimeDistortionDummy=startingTimeDistortion+round(widthSHField*cosd(180-orientation)/velocityDistortion);
end;

% loop for number of evaporator loops
for k=1:nString

startingTimeDistortion=startingTimeDistortionDummy+round((k-1)*widthEvapField/(nString-1)*sind(270-orientation)/velocityDistortion);

dims1=discIrr*2+2;
dims2=discIrr+1;


checkLengthDistortion=ceil((lengthDistortion-lengthTube/2)/(lengthTube/discIrr));
discDistortion=ceil(lengthDistortion/lengthTube*discIrr);
%
%
%

if (orientation >= 180)
%  orientation West
if checkLengthDistortion < 0
    for j=1:nDistortion

    for i=2:(dims1-2)/4+1

     dummystring=round(startingTimeDistortion+(i-2)*lengthTubeHor/velocityDistortion+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
    dummystring1=round(startingTimeDistortion+(i-2)*lengthTubeHor/velocityDistortion+deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
    S2((i-2)*2+2+(j-1)*((((dims1-2)/4+1)-2)*2+2),1)={num2str(dummystring)};
    S2((i-2)*2+3+(j-1)*((((dims1-2)/4+1)-2)*2+2),1)={num2str(dummystring1)};
    end
    end

    for j=1:nDistortion
         for i=1:(dims1-6)/4+1
     dummystring=round(startingTimeDistortion+(i-1)*lengthTubeHor/velocityDistortion+deltaTimeDistortion+deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
     dummystring1=round(startingTimeDistortion+(i-1)*lengthTubeHor/velocityDistortion+deltaTimeDistortion+2*deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
     S2((i-1)*2+2+nDistortion*(discIrr)+(j-1)*discIrr,1)={num2str(dummystring)};
    S2((i-1)*2+3+nDistortion*discIrr+(j-1)*discIrr,1)={num2str(dummystring1)};
         end
    end

    S2(1,1)={'0'};
    S2(end+1,1)={num2str(endTimeSim)};
    S2=unique(S2(:,1));
    S2(1,1)={'0'};
    S2(1:end,2:dims2)={num2str(Irr)};


    for h=1:nDistortion
    for i=2:discIrr/2+1
        for j=1:length(S2(:,1))
        if (str2num(S2{j,1})>=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion))) && (str2num(S2{j,1})<=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+deltaTimeDistortion+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion)));
            S2(j,discIrr/2-i+3)={num2str(Irr-deltaIrr)};
        end
        end
         for j=1:length(S2(:,1))
        if (str2num(S2{j,1})>=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion))) && (str2num(S2{j,1})<=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+deltaTimeDistortion+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion)));
            S2(j,discIrr/2+i)={num2str(Irr-deltaIrr)};
        end
         end
    end
    end
else
    for j=1:nDistortion

    for i=2:(dims1-2)/4+1

     dummystring=round(startingTimeDistortion+(i-2)*lengthTubeHor/velocityDistortion+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
    dummystring1=round(startingTimeDistortion+(i-2)*lengthTubeHor/velocityDistortion+deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
    S2((i-2)*2+2+(j-1)*((((dims1-2)/4+1)-2)*2+2),1)={num2str(dummystring)};
    S2((i-2)*2+3+(j-1)*((((dims1-2)/4+1)-2)*2+2),1)={num2str(dummystring1)};
    end
    end

    for j=1:nDistortion
         for i=1:(dims1-6)/4+1
     dummystring=round(startingTimeDistortion+(i-1)*lengthTubeHor/velocityDistortion+deltaTimeDistortion+deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
     dummystring1=round(startingTimeDistortion+(i-1)*lengthTubeHor/velocityDistortion+deltaTimeDistortion+2*deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
     S2((i-1)*2+2+nDistortion*(discIrr)+(j-1)*discIrr,1)={num2str(dummystring)};
    S2((i-1)*2+3+nDistortion*discIrr+(j-1)*discIrr,1)={num2str(dummystring1)};
         end
    end

         S2(1,1)={'0'};
    S2(end+1,1)={num2str(endTimeSim)};
    S2=unique(S2(:,1));
    S2(1,1)={'0'};
    S2(1:end,2:dims2)={num2str(Irr)};


      for h=1:nDistortion
    for i=2:discIrr/2+1
        for j=1:length(S2(:,1))
        if (str2num(S2{j,1})>=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion))) && (str2num(S2{j,1})<=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+deltaTimeDistortion+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion)));
            S2(j,discIrr/2-i+3)={num2str(Irr-deltaIrr)};
        end
        end
         for j=1:length(S2(:,1))
        if (str2num(S2{j,1})>=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion))) && (str2num(S2{j,1})<=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+deltaTimeDistortion+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion)));
            S2(j,discIrr/2+i)={num2str(Irr-deltaIrr)};
        end
        end
    end
      end

end

%   orientation East
else
    for j=1:nDistortion

    for i=2:(dims1-2)/4+1

     dummystring=round(startingTimeDistortion+(i-2)*lengthTubeHor/velocityDistortion+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
    dummystring1=round(startingTimeDistortion+(i-2)*lengthTubeHor/velocityDistortion+deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
    S2((i-2)*2+2+(j-1)*((((dims1-2)/4+1)-2)*2+2),1)={num2str(dummystring)};
    S2((i-2)*2+3+(j-1)*((((dims1-2)/4+1)-2)*2+2),1)={num2str(dummystring1)};
    end
    end

    for j=1:nDistortion
         for i=1:(dims1-6)/4+1
     dummystring=round(startingTimeDistortion+(i-1)*lengthTubeHor/velocityDistortion+deltaTimeDistortion+deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
     dummystring1=round(startingTimeDistortion+(i-1)*lengthTubeHor/velocityDistortion+deltaTimeDistortion+2*deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
     S2((i-1)*2+2+nDistortion*(discIrr)+(j-1)*discIrr,1)={num2str(dummystring)};
    S2((i-1)*2+3+nDistortion*discIrr+(j-1)*discIrr,1)={num2str(dummystring1)};
         end
    end

    S2(1,1)={'0'};
    S2(end+1,1)={num2str(endTimeSim)};
    S2=unique(S2(:,1));
    S2(1,1)={'0'};
    S2(1:end,2:dims2)={num2str(Irr)};


    for h=1:nDistortion
        for i=2:discIrr/2+1
        for j=1:length(S2(:,1))
        if (str2num(S2{j,1})>=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion))) && (str2num(S2{j,1})<=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+deltaTimeDistortion+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion)));
            S2(j,i)={num2str(Irr-deltaIrr)};
        end
        end
         for j=1:length(S2(:,1))
        if (str2num(S2{j,1})>=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion))) && (str2num(S2{j,1})<=round(startingTimeDistortion+((i-2)*lengthTubeHor)/velocityDistortion+deltaTime+deltaTimeDistortion+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion)));
            S2(j,discIrr-i+3)={num2str(Irr-deltaIrr)};
        end
        end
    end
    end
% end

end


%% create new txt file

    name =char(strcat('EvaporatorLoop', num2str(k)))   ;            % name of file
    pthstr ='C:/Work/Clara_Solar_InputData/InputIrradiation';    % path of file

gile=[pthstr '/' name '.txt']; ... Modelica Input File

%%%%%%%%%Modelica Input File%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen(gile, 'wt');
fprintf(fid,'%s\n','#1');
%fprintf(fid,'%s\n',['#Begin of Data at ' R{1}]);
%fprintf(fid,'%s\n',['#End of Data at ' R{2}]);
fprintf(fid,'%s\n',['double dat(' num2str(length(S2(:,1))) ',' num2str(length(S2(1,:))) ')']);
cnt=1:dims2;
fprintf(fid,'%s\n',['# ' num2str(cnt)]);
%fprintf(fid,'%s\n',['#' strhcat(strrep(N2,' ','_'))]);
%fprintf(fid,'%s\n',['#' strhcat(strrep(T,' ','_'))]);
for i=1:length(S2(:,1))
    S2N = strjoin(strcat(S2(i,:)), ' ');
    %fprintf(fid,'%s \n',strcat(S2{i,:}));
    fprintf(fid,'%s\n',S2N);
end
fprintf(fid,'%s\n',['#Collector number' num2str(k)]);
fprintf(fid,'%s\n',['#Date of creation:' date]);
fprintf(fid,'%s\n',['#lengthTube=' num2str(lengthTube) ' discIrr=' num2str(discIrr) ' lengthDistortion=' num2str(lengthDistortion) ' velocityDistortion=' num2str(velocityDistortion) ' endTimeSim=' num2str(endTimeSim) ' deltaTimeDistortion=' num2str(deltaTimeDistortion) ' deltaTime=' num2str(deltaTime) ' Irr=' num2str(Irr) ' deltaIrr=' num2str(deltaIrr) ' startingTimeDistortion=' num2str(startingTimeDistortion) ]);

fclose(fid);

end;

%% create table data superheater

if (orientation <= 90) || (orientation >= 270)
    startingTimeDistortionDummy=startingTimeDistortionDummy+round(widthSHField*sind(90-orientation)/velocityDistortion);
else
    startingTimeDistortionDummy=startingTimeDistortionDummy2;
end;


% loop for number of superheater loops
for k=1:nStringSH

    startingTimeDistortion=startingTimeDistortionDummy+round((k-1)*widthSHField/(nStringSH-1)*sind(270-orientation)/velocityDistortion);

dims1=discIrrSH*4+2;
dims2=discIrrSH+1;


checkLengthDistortion=ceil((lengthDistortion-lengthTubeSH/2)/(lengthTubeSH/discIrrSH));
discDistortion=ceil(lengthDistortion/lengthTubeSH*discIrrSH);
%
%
%
if (orientation >= 180)
%  orientation West


% if checkLengthDistortion < 0
%     for i=2:(dims1-2)/2+1
%      dummystring=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion);
%     dummystring1=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion+deltaTime);
%     S2((i-2)*2+2,1)={num2str(dummystring)};
%     S2((i-2)*2+3,1)={num2str(dummystring1)};
%
%     end
%          for i=(dims1-2)/2+2:(dims1)/2
%      dummystring=round(startingTimeDistortion+((i-((dims1-2)/4+2))*lengthTubeHorSH)/velocityDistortion+deltaTimeDistortion+deltaTime);
%      dummystring1=round(startingTimeDistortion+((i-((dims1-2)/4+2))*lengthTubeHorSH)/velocityDistortion+deltaTimeDistortion+2*deltaTime);
%     S2((i-2)*2+2,1)={num2str(dummystring)};
%     S2((i-2)*2+3,1)={num2str(dummystring1)};
%          end
%
%     S2(1,1)={'0'};
%     S2(end+1,1)={num2str(endTimeSim)};
%     S2=unique(S2(:,1));
%     S2(1,1)={'0'};
%     S2(1:end,2:dims2)={num2str(Irr)};
%
%
%     for i=2:discIrrSH+1
%         for j=1:length(S2(:,1))
%         if (str2num(S2{j,1})>=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion+deltaTime)) && (str2num(S2{j,1})<=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion+deltaTime+deltaTimeDistortion))
%             S2(j,discIrr/2-i+3)={num2str(Irr-deltaIrr)};
%         end
%         end
%     end
%
% else
        for j=1:nDistortion

    for i=2:(dims1-2)/4+1

     dummystring=round(startingTimeDistortion+(i-2)*lengthTubeHorSH/velocityDistortion+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
    dummystring1=round(startingTimeDistortion+(i-2)*lengthTubeHorSH/velocityDistortion+deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
    S2((i-2)*2+2+(j-1)*((((dims1-2)/4+1)-2)*2+2),1)={num2str(dummystring)};
    S2((i-2)*2+3+(j-1)*((((dims1-2)/4+1)-2)*2+2),1)={num2str(dummystring1)};
    end
    end

    for j=1:nDistortion
         for i=1:(dims1-6)/4+1
     dummystring=round(startingTimeDistortion+(i-1)*lengthTubeHorSH/velocityDistortion+deltaTimeDistortion+deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
     dummystring1=round(startingTimeDistortion+(i-1)*lengthTubeHorSH/velocityDistortion+deltaTimeDistortion+2*deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion));
     S2((i-1)*2+2+nDistortion*(discIrrSH*2)+(j-1)*discIrrSH*2,1)={num2str(dummystring)};
    S2((i-1)*2+3+nDistortion*discIrrSH*2+(j-1)*discIrrSH*2,1)={num2str(dummystring1)};
         end
    end

    S2(1,1)={'0'};
    S2(end+1,1)={num2str(endTimeSim)};
    S2=unique(S2(:,1));
    S2(1,1)={'0'};
    S2(1:end,2:dims2)={num2str(Irr)};


    for h=1:nDistortion
    for i=2:discIrrSH+1
        for j=1:length(S2(:,1))
        if (str2num(S2{j,1})>=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion+deltaTime+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion))) && (str2num(S2{j,1})<=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion+deltaTime+deltaTimeDistortion)+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion))
            S2(j,i)={num2str(Irr-deltaIrr)};
        end
        end
    end
    end

% end

%   orientation East
else
% if checkLengthDistortion < 0
%     for i=2:(dims1-2)/2+1
%      dummystring=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion+lengthTubeHorSHPowerBlock/velocityDistortion);
%     dummystring1=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion+lengthTubeHorSHPowerBlock/velocityDistortion+deltaTime);
%     S2((i-2)*2+2,1)={num2str(dummystring)};
%     S2((i-2)*2+3,1)={num2str(dummystring1)};
%
%     end
%          for i=(dims1-2)/2+2:(dims1)/2
%      dummystring=round(startingTimeDistortion+((i-((dims1-2)/4+2))*lengthTubeHorSH)/velocityDistortion+deltaTimeDistortion+lengthTubeHorSHPowerBlock/velocityDistortion+deltaTime);
%      dummystring1=round(startingTimeDistortion+((i-((dims1-2)/4+2))*lengthTubeHorSH)/velocityDistortion+deltaTimeDistortion+lengthTubeHorSHPowerBlock/velocityDistortion+2*deltaTime);
%     S2((i-2)*2+2,1)={num2str(dummystring)};
%     S2((i-2)*2+3,1)={num2str(dummystring1)};
%          end
%
%     S2(1,1)={'0'};
%     S2(end+1,1)={num2str(endTimeSim)};
%     S2=unique(S2(:,1));
%     S2(1,1)={'0'};
%     S2(1:end,2:dims2)={num2str(Irr)};
%
%     for i=2:discIrrSH+1
%         for j=1:length(S2(:,1))
%         if (str2num(S2{j,1})>=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion+lengthTubeHorSHPowerBlock/velocityDistortion+deltaTime)) && (str2num(S2{j,1})<=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion+lengthTubeHorSHPowerBlock/velocityDistortion+deltaTime+deltaTimeDistortion))
%             S2(j,discIrrSH+3-i)={num2str(Irr-deltaIrr)};
%         end
%         end
%     end
%
% else
%
        for j=1:nDistortion

    for i=2:(dims1-2)/4+1

     dummystring=round(startingTimeDistortion+(i-2)*lengthTubeHorSH/velocityDistortion+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion)+lengthTubeHorSHPowerBlock/velocityDistortion);
    dummystring1=round(startingTimeDistortion+(i-2)*lengthTubeHorSH/velocityDistortion+deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion)+lengthTubeHorSHPowerBlock/velocityDistortion);
    S2((i-2)*2+2+(j-1)*((((dims1-2)/4+1)-2)*2+2),1)={num2str(dummystring)};
    S2((i-2)*2+3+(j-1)*((((dims1-2)/4+1)-2)*2+2),1)={num2str(dummystring1)};
    end
    end

    for j=1:nDistortion
         for i=1:(dims1-6)/4+1
     dummystring=round(startingTimeDistortion+(i-1)*lengthTubeHorSH/velocityDistortion+deltaTimeDistortion+deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion)+lengthTubeHorSHPowerBlock/velocityDistortion);
     dummystring1=round(startingTimeDistortion+(i-1)*lengthTubeHorSH/velocityDistortion+deltaTimeDistortion+2*deltaTime+(j-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion)+lengthTubeHorSHPowerBlock/velocityDistortion);
     S2((i-1)*2+2+nDistortion*(discIrrSH*2)+(j-1)*discIrrSH*2,1)={num2str(dummystring)};
    S2((i-1)*2+3+nDistortion*discIrrSH*2+(j-1)*discIrrSH*2,1)={num2str(dummystring1)};
         end
    end

    S2(1,1)={'0'};
    S2(end+1,1)={num2str(endTimeSim)};
    S2=unique(S2(:,1));
    S2(1,1)={'0'};
    S2(1:end,2:dims2)={num2str(Irr)};


        for h=1:nDistortion
        for i=2:discIrrSH+1
        for j=1:length(S2(:,1))
        if (str2num(S2{j,1})>=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion+lengthTubeHorSHPowerBlock/velocityDistortion+deltaTime+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion))) && (str2num(S2{j,1})<=round(startingTimeDistortion+((i-2)*lengthTubeHorSH)/velocityDistortion+lengthTubeHorSHPowerBlock/velocityDistortion+deltaTime+deltaTimeDistortion+(h-1)*(diffLengthDistortion/velocityDistortion+deltaTimeDistortion)));
            S2(j,discIrrSH+3-i)={num2str(Irr-deltaIrr)};
        end
        end
    end
        end


end


%% create new txt file

    name =char(strcat('SuperheaterLoop', num2str(k)));               % name of file
    pthstr ='C:/Work/Clara_Solar_InputData/InputIrradiation';    % path of file

gile=[pthstr '/' name '.txt']; ... Modelica Input File

%%%%%%%%%Modelica Input File%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen(gile, 'wt');
fprintf(fid,'%s\n','#1');
%fprintf(fid,'%s\n',['#Begin of Data at ' R{1}]);
%fprintf(fid,'%s\n',['#End of Data at ' R{2}]);
fprintf(fid,'%s\n',['double dat(' num2str(length(S2(:,1))) ',' num2str(length(S2(1,:))) ')']);
cnt=1:dims2;
fprintf(fid,'%s\n',['# ' num2str(cnt)]);
%fprintf(fid,'%s\n',['#' strhcat(strrep(N2,' ','_'))]);
%fprintf(fid,'%s\n',['#' strhcat(strrep(T,' ','_'))]);
for i=1:length(S2(:,1))
    S2N = strjoin(strcat(S2(i,:)), ' ');
    %fprintf(fid,'%s \n',strcat(S2{i,:}));
    fprintf(fid,'%s\n',S2N);
end
fprintf(fid,'%s\n',['#Collector number' num2str(k)]);
fprintf(fid,'%s\n',['#Date of creation:' date]);
fprintf(fid,'%s\n',['#lengthTubeSH=' num2str(lengthTubeSH) ' discIrrSH=' num2str(discIrrSH) ' lengthDistortion=' num2str(lengthDistortion) ' velocityDistortion=' num2str(velocityDistortion) ' endTimeSim=' num2str(endTimeSim) ' deltaTimeDistortion=' num2str(deltaTimeDistortion) ' deltaTime=' num2str(deltaTime) ' Irr=' num2str(Irr) ' deltaIrr=' num2str(deltaIrr) ' startingTimeDistortion=' num2str(startingTimeDistortion) ]);

fclose(fid);

end;
