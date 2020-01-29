function haxesL=LogoIEO(Type,Position)

if nargin==0
    Type=1;
    Position=[0.58 0.38 0.14 0.14];
elseif nargin==1
    switch Type
        case 1
            Position=[0.58 0.35 0.14 0.14];
        case 2
            Position=[0.58 0.35 0.20 0.20];
        case 3
            Position=[0.55 0.08 0.25 .14];
        case 4
            Position=[0.55 0.08 0.25 .14];
        case 5
            Position=[0.58 0.35 0.14 0.14];
    end
end

LogoFile=strcat('LogoIEO',num2str(Type),'.png');

[img, map, alphachannel]=imread(LogoFile);

haxesL=axes;
haxesL.PlotBoxAspectRatio=[1 1 1];

himagL=image(img, 'AlphaData', alphachannel);

haxesL.Position=Position;
haxesL.Color='none';
haxesL.XTickLabel='';
haxesL.YTickLabel='';
haxesL.Box='off';
haxesL.DataAspectRatioMode='auto';
haxesL.DataAspectRatio=[1 1 1];
haxesL.XColor='none';
haxesL.YColor='none';

end

