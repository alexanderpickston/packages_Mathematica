(* ::Package:: *)

BeginPackage["BeamProfiling`"]


Remove[BeamProfiling,CheckBeamOverlap,CheckBeamInFocus];


BeamProfiling::usage = 
"Visualize, track and use pointer on beam.
The argument is the RasterSize of the Capture: {324,243} is the default value.";

CheckBeamOverlap::usage=
"Check whether different beams are overlapped.
The first argument is the threshold of the image binarize: it's default value is 0.2.
CheckBeamOverlap takes as second (optional) argument a list with the file names.
It prints several data 
If called with no second argument CheckBeamOverlap checks all the jpgs and pngs files in the notebook directory.
"

CheckBeamInFocus::usage=
"Check the dimension and the centre position of the beam focused.
The first argument is an approximation of the beam width, used to crop the data. The default value is 200 and it usually works for focused beams..
If called with no second argument CheckBeamInFocus checks all the jpgs and pngs files in the notebook directory.
As input pictures, use the 2592H x 1944V images from the CMOS camera with 2.2\[Mu]m pixel.
";


Begin["`Private`"]


BeamProfiling[raster_:{324,243}]:=

Block[{pts,target,edge,tracking,img,camera,fit,dimension,fileName,dirName},
fitList1=ConstantArray[None,200];
fitList2=ConstantArray[None,200];
Manipulate[
LocatorPane[
Dynamic[
pts,

With[{update=#[[1]]-pts[[1]]},
If[Max@Abs[update]==0,
pts=#,
pts+={update,update}]]&],

Dynamic@Show[
img=CurrentImage[RasterSize->raster,ImagingDevice->camera],

If[target==True,
Graphics[{Blue,Thick,Circle[pts[[1]],Abs[Subtract@@pts]*Sqrt[2]]}]
,Unevaluated[Sequence[]]
]
,
If[tracking==True,
edge=EdgeDetect[Binarize[img,0.2]];
fit=circfit[N@PixelValuePositions[edge,1]][[1]];
AppendTo[fitList1,fit[[2,2,1]]];
AppendTo[fitList2,fit[[2,2,2]]];
Drop[fitList1,1];
Drop[fitList2,1];

Show[
ColorReplace[edge,{White->Red,Black->GrayLevel[0, 0]}]
,Graphics[{Red,Thick,Circle[fit[[2,2]],fit[[3,2]]]
,Inset[Style[TextGrid[{{"Centre: "<>ToString@Round@fit[[2,2]] },{"Radius: "<>ToString@NumberForm[Abs@fit[[3,2]] ,{3,1}]}}],20,Red],{70,220}]
,Inset[ListPlot[{fitList1[[-200;;]]},ImageSize->130,PlotStyle->LightGreen,PlotRange->All,Frame->True,FrameStyle->White,LabelStyle->White],{280,210}]
,Inset[ListPlot[{fitList2[[-200;;]]},ImageSize->130,PlotStyle->LightGreen,PlotRange->All,Frame->True,FrameStyle->White,LabelStyle->White],{280,160}]
}
]
,Graphics[{Red,PointSize[0.01],Point[fit[[2,2]]]}]
]
,Unevaluated[Sequence[]]
]

,ImageSize->dimension
]

,If[target==True,Appearance->{Style["\[CircleDot]",Blue,30],Style["\[FilledCircle]",Blue,16]},Appearance->{Style[Black,Opacity[0]],Style[Black,Opacity[0]]}]

]
,{{pts,{raster/2.,raster/1.7}},ControlType->None}

,Pane[Grid[
{
{SetterBar[Dynamic@camera,$ImagingDevices],SpanFromLeft},
{Control@{{dimension,600,"Image Size"},20,1000,ControlType->Slider},SpanFromLeft},
{Control@{{target,False,"Show Target"},{True,False}},Button[Style[Text["Reset Target"],14],pts={{324,243}/2.,{281,200}/2.},ImageSize->{Automatic,30}]},
{Control@{{tracking,False,"Show Tracking"},{True,False}},Button[Style[Text["Target On Tracking"],14],pts={fit[[2,2]],fit[[2,2]]+{fit[[3,2]]/Sqrt[2],fit[[3,2]]/Sqrt[2]}},ImageSize->{Automatic,30}]},
{Control@{{fileName,"img","Set File Name"},ControlType->InputField[ToString@HoldForm@#&,String]},SpanFromLeft},
{Control@{{dirName,NotebookDirectory[],"Set Directory"},ControlType->InputField[ToString@HoldForm@#&,String]},SpanFromLeft},{Button[Style[Text["Snapshot"],14],Export[ToString@dirName<>ToString@fileName<>".png",img],ImageSize->{Automatic,30}],SpanFromLeft}
}
]]
,LabelStyle->14
]

]


CheckBeamOverlap[threshold_:0.2,files__:Hold[FileNames[NotebookDirectory[]<>"*.png"]~Join~FileNames[NotebookDirectory[]<>"*.jpg"]]]:=
Block[{pics,filePaths,edges,fits,colors,names,l,binar,dims},
filePaths=Flatten[{ReleaseHold@files}];

names=FileBaseName/@filePaths;
l=Length[filePaths];
pics=Import[#]&/@filePaths;
binar=Binarize[#,threshold]&/@pics;
dims=ImageData[binar[[1]]]//Dimensions;

edges=EdgeDetect/@binar;
fits=circfit[N@PixelValuePositions[#,1]][[1]]&/@edges;

If[l<=10,
colors={Red,Blue,Green,Magenta,Orange,Purple,Yellow,Pink,Cyan,Brown},
colors=RandomColor[l]
];

Print[
Show[
Graphics[Style[{colors[[#]],Point["center"]/.fits[[#]]}]&/@Range[l]],
Graphics[Style[{colors[[#]],Circle["center","radius"]/.fits[[#]]}]&/@Range[l]]
,ImageSize->Medium]
];

Print[
TableForm[Style[TextGrid[
{{names[[#]],NumberForm["center",{4,1}],NumberForm["radius",{4,1}]}}/.fits[[#]]],14,colors[[#]]]&/@Range[l]]];

Print[
GraphicsGrid[
{Flatten[
Show[
pics[[#]],
Graphics[{Style[{colors[[#]],Circle["center","radius"]/.fits[[#]]}]
,Inset[Style[TextGrid[{{names[[#]]}}],14,colors[[#]]],dims*{0.33,0.66}]}
]]&/@Range[l]
]
}
,ImageSize->1000]
];

{names[[#]],
"center",
"radius",
Show[
pics[[#]],
Graphics[{Style[{colors[[#]],Circle["center","radius"]/.fits[[#]]}]
,Inset[Style[TextGrid[{{names[[#]]}}],14,colors[[#]]],dims*{0.33,0.66}]}
]]}/.fits[[#]]&/@Range[l]

]


CheckBeamInFocus[widthbeam_:200,files__:Hold[FileNames[NotebookDirectory[]<>"*.png"]~Join~FileNames[NotebookDirectory[]<>"*.jpg"]]]:=
Block[{a2D1,x2D1,y2D1,\[Sigma]2D1H,\[Sigma]2D1V,pics,filePaths,pixelsize,names,l,max,cutImg,maxCoord,fit2D,center,FWHMHorizontal2D,FWHMVertical2D},
filePaths=Flatten[{ReleaseHold@files}];

names=FileBaseName/@filePaths;
l=Length[filePaths];

pixelsize=2.2;

pics=Import[#,"GrayLevels"]&/@filePaths;
max=Mean[Position[#,Max[#]]]&/@pics;

cutImg=
Flatten[
Table[
{j*pixelsize,i*pixelsize,Total[pics[[#,i;;i,j;;j]],2]}
,{i,Round[max[[#,1]]-widthbeam],Round[max[[#,1]]+widthbeam],1}
,{j,Round[max[[#,2]]-widthbeam],Round[max[[#,2]]+widthbeam],1}
],1]&/@Range[l];

cutImg//Dimensions;

Print[GraphicsGrid[{ListDensityPlot[#[[;;;;10]],PlotRange->All,ImageSize->Small]&/@cutImg}]];

maxCoord=SortBy[#,Last][[-1]]&/@cutImg;

fit2D=
FindFit[cutImg[[#]],a2D1*Exp[-((x-x2D1)^2/(2\[Sigma]2D1H^2))]*Exp[-((y-y2D1)^2/(2\[Sigma]2D1V^2))]
,{{a2D1,maxCoord[[#,3]]},{x2D1,maxCoord[[#,1]]},{y2D1,maxCoord[[#,2]]},{\[Sigma]2D1H,100},{\[Sigma]2D1V,100}},{x,y}]&/@Range[l];

center={x2D1,y2D1}/.fit2D;
FWHMHorizontal2D=2Sqrt[2Log[2]]*Abs@\[Sigma]2D1H/.fit2D;
FWHMVertical2D=2Sqrt[2Log[2]]*Abs@\[Sigma]2D1V/.fit2D;

Print[
TableForm[
Prepend[Transpose@{names,center,FWHMHorizontal2D,
FWHMVertical2D,(FWHMHorizontal2D+FWHMVertical2D)/2,
(FWHMHorizontal2D+FWHMVertical2D)/Sqrt[2Log[2]]},{"Path","(\!\(\*SubscriptBox[\(x\), \(0\)]\),\!\(\*SubscriptBox[\(y\), \(0\)]\))","FWHM_X","FWHM_Y","Avarage","waist (\!\(\*FractionBox[\(1\), SuperscriptBox[\(\[ExponentialE]\), \(2\)]]\))"}]
]
];

Transpose@{names,center,FWHMHorizontal2D,
FWHMVertical2D,(FWHMHorizontal2D+
FWHMVertical2D)/2}

]


circfit[pts_]:=Block[{reg,lm,bf,exp,center,rad,x,y},
reg={2 #1,2 #2,#2^2+#1^2}&@@@pts;
lm=LinearModelFit[reg,{1,x,y},{x,y}];
bf=lm["BestFitParameters"];
exp=(x-#2)^2+(y-#3)^2-#1-#2^2-#3^2&@@bf;
{center,rad}={{#2,#3},Sqrt[#2^2+#3^2+#1]}&@@bf;
circlefit[{"expression"->exp,"center"->center,"radius"->rad}]];circlefit[list_][field_]:=field/.list;
circlefit[list_]["Properties"]:=list/.Rule[field_,_]:>field;
circlefit/:ReplaceAll[fields_,circlefit[list_]]:=fields/.list;
Format[circlefit[list_],StandardForm]:=HoldForm[circlefit]["<"<>ToString@Length@list<>">"]


End[]


EndPackage[]
