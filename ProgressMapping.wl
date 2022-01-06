(* ::Package:: *)

(* ::Title:: *)
(*Progress Mapping*)


(* ::Chapter:: *)
(*Begin package*)


BeginPackage["ProgressMapping`"]


(* ::Section::Closed:: *)
(*Package description*)


ProgressMap::usage =
 "ProgressMap[f,expr] is a Map implementation with progress bar. levelspec is always {1}.
\"ShowInfo\"\[Rule]True show the detailed version of the progress. \"Label\"->\"Custom label\" shows a custom description label.";

ProgressParallelMap::usage =
 "ProgressParallelMap[f,expr] is a ParallelMap implementation with progress bar. levelspec is always {1}. 
\"ShowInfo\"\[Rule]True show the detailed version of the progress bar. \"Label\"->\"Custom label\" shows a custom description label. It inherits all options of ParallelMap.";

ProgressMapThread::usage = 
"ProgressMapThread[f, {{\!\(\*SubscriptBox[
StyleBox[\"a\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"a\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\),\!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\)},{\!\(\*SubscriptBox[
StyleBox[\"b\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"b\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\),\!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\)},\!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\)}] is a MapThread implementation with progress bar. levelspec is always {1}. 
\"ShowInfo\"\[Rule]True show the detailed version of the progress bar. \"Label\"->\"Custom label\" shows a custom description label."

ProgressParallelMapThread::usage = 
"ProgressParallelMapThread[f, {{\!\(\*SubscriptBox[
StyleBox[\"a\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"a\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\),\!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\)},{\!\(\*SubscriptBox[
StyleBox[\"b\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"b\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\),\!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\)},\!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\)}] is a ParallelMapThread implementation with progress bar. levelspec is always {1}. 
\"ShowInfo\"\[Rule]True show the detailed version of the progress bar. \"Label\"->\"Custom label\" shows a custom description label."
  
ProgressTable::usage =
 "Table implementation with progress bar. Same usage as Table. 
\"ShowInfo\"\[Rule]True show the detailed version of the progress bar. \"Label\"->\"Custom label\" shows a custom description label.";

ProgressParallelTable::usage =
"ParallelTable implementation with progress bar. Same usage as ParallelTable. 
\"ShowInfo\"\[Rule]True show the detailed version of the progress bar. \"Label\"->\"Custom label\" shows a custom description label.
It inherits all the options from ParallelTable.";

ProgressNest::usage =
 "Nest implementation with progress bar. Same usage as Nest. 
\"ShowInfo\"\[Rule]True show the detailed version of the progress bar. \"Label\"->\"Custom label\" shows a custom description label.";

ProgressNestList::usage =
 "NestList implementation with progress bar. Same usage as NestList. 
\"ShowInfo\"\[Rule]True show the detailed version of the progress bar. \"Label\"->\"Custom label\" shows a custom description label.";


DefaultIndicator::usage = "DefaultIndicator[indexProgress, totalSize] show the default indicator.";
DetailedIndicator::usage = "DetailedIndicator[indexProgress, totalSize, startTime, label] show the detailed indicator. Can be used in custom Monitors.";


(* ::Section::Closed:: *)
(*Error messages*)


ProgressParallelMap::exprlengtherr = "expr must have at least a length greater than 1";
ProgressMap::exprlengtherr = "expr must have at least a length greater than 1";
ProgressMapThread::exprlengtherr = "expr must have at least a length greater than 1";
ProgressMapThread::exprdeptherr = "expr must have a depth of at least 3, current depth is `1`";
ProgressParallelMapThread::exprlengtherr = "expr must have at least a length greater than 1";
ProgressParallelMapThread::exprdeptherr = "expr must have a depth of at least 3, current depth is `1`";


(* ::Chapter:: *)
(*Begin implementation*)


Begin["`Private`"]


(* ::Section::Closed:: *)
(*Clock*)


(* ::Text:: *)
(*Simple hh/mm/ss implementation of a clock.*)


GetSeconds[time_] := IntegerString[Round[Mod[time, 60]], 10, 2];
GetMinutes[time_]:= IntegerString[Mod[Floor[time/60], 60], 10, 2];
GetHours[time_] := IntegerString[Floor[time/3600], 10, 2];
ClockFormat[time_] := StringJoin[GetHours[time], ":", GetMinutes[time], ":", GetSeconds[time]];


(* ::Section::Closed:: *)
(*Progress indicator*)


(* ::Text:: *)
(*Progress indicator to be used by ProgressMap and ProgressTable family of functions.*)


DefaultIndicator[indexProgress_, totalSize_] := ProgressIndicator[indexProgress, {1, totalSize}];

DetailedIndicator[indexProgress_, totalSize_, startTime_, label_:"Evaluating..."] := 
Block[{progressString, remainingTime, remainingTimeString, indicator, ellapsedTimeString},
	progressString = Row[{Style["Progress: ", Bold], ToString[indexProgress], "/", ToString[totalSize]}];
	ellapsedTimeString = Row[{Style["Elapsed time: ", Bold], ClockFormat[AbsoluteTime[] - startTime]}];

	If[indexProgress != 0,
		remainingTime = ((AbsoluteTime[] - startTime) / indexProgress)*(totalSize - indexProgress);
		remainingTimeString = Row[{Style["Remaining: ", Bold], ClockFormat[remainingTime]}];
		,
		remainingTimeString = Row[{Style["Remaining: ", Bold], "Unknown"}];
	];

	indicator = Panel[
		Column[
			{
				Style[label, Bold],
				DefaultIndicator[indexProgress, totalSize],
				progressString,
				ellapsedTimeString,
				remainingTimeString
			}
		]
	];

	Return[indicator];
];


(* ::Section::Closed:: *)
(*Maps with progress*)


SetAttributes[ProgressParallelMap, HoldFirst];
ProgressParallelMap[f_, expr_, opts: OptionsPattern[{"ShowInfo"->True, "Label"->"Evaluating...", Parallelize}]] :=
Module[{startTime = AbsoluteTime[], indexProgress = 0, output},
	SetSharedVariable[indexProgress];
	If[Length[expr]<1, Message[ProgressParallelMap::exprlengtherr];Return[$Failed]];

	Monitor[
		ParallelMap[
			(
				output = f[#];
				indexProgress++;
				output
			)&,
			expr,
			FilterRules[{opts}, Options[Parallelize]]
		]
		,
		If[OptionValue["ShowInfo"],
			DetailedIndicator[indexProgress, Length[expr], startTime, OptionValue["Label"]]
			,
			DefaultIndicator[indexProgress, Length[expr]]
		]
	]
]; 

SetAttributes[ProgressMap, HoldFirst];
ProgressMap[f_, expr_, OptionsPattern[{"ShowInfo"->True, "Label"->"Evaluating..."}]] :=
Module[{startTime = AbsoluteTime[], indexProgress = 0, output},
	Monitor[
		If[Length[expr]<1, Message[ProgressMap::exprlengtherr];Return[$Failed]];
		Map[
			(
				output = f[#];
				indexProgress++;
				output
			)&,
			expr
		]
		,
		If[OptionValue["ShowInfo"],
			DetailedIndicator[indexProgress, Length[expr], startTime, OptionValue["Label"]]
			,
			DefaultIndicator[indexProgress, Length[expr]]
		]
	]
];


(* ::Section::Closed:: *)
(*MapThread with progress*)


SetAttributes[ProgressMapThread,HoldFirst];
ProgressMapThread[f_, expr_, OptionsPattern[{"ShowInfo"->True,"Label"->"Evaluating..."}]] :=
Module[{startTime = AbsoluteTime[], indexProgress = 0, output},
	If[Length[expr] < 1, Message[ProgressMapThread::exprlengtherr];Return[$Failed]];
	If[Depth[expr] < 3, Message[ProgressMapThread::exprdeptherr, Depth[expr]];Return[$Failed]];

	Monitor[
		MapThread[
			(
				output = f[##];
				indexProgress++;
				output
			)&
			,
			expr
		]
		,
		If[OptionValue["ShowInfo"],
			DetailedIndicator[indexProgress, Length[First[expr]], startTime, OptionValue["Label"]]
			,
			DefaultIndicator[indexProgress, Length[First[expr]]]
		]
	]
]; 

SetAttributes[ProgressParallelMapThread,HoldFirst];
ProgressParallelMapThread[f_, expr_, OptionsPattern[{"ShowInfo"->True,"Label"->"Evaluating..."}]] :=
Module[{startTime = AbsoluteTime[], indexProgress = 0, output},
	SetSharedVariable[indexProgress];
	If[Length[expr] < 1, Message[ProgressParallelMapThread::exprlengtherr];Return[$Failed]];
	If[Depth[expr] < 3, Message[ProgressParallelMapThread::exprdeptherr, Depth[expr]];Return[$Failed]];

	Monitor[
		ParallelMapThread[
		(
			output = f[##];
			indexProgress++;
			output
		)&
		,
		expr
		]
		,
		If[OptionValue["ShowInfo"],
			DetailedIndicator[indexProgress, Length[First[expr]], startTime, OptionValue["Label"]]
			,
			DefaultIndicator[indexProgress, Length[First[expr]]]
		]
	]
]; 


(* ::Section::Closed:: *)
(*Tables with progress*)


SetAttributes[ProgressTable, HoldFirst];
Options[ProgressTable] = {"ShowInfo"->True, "Label"->"Evaluating..."};

ProgressTable[expr_, n_Integer, opts: OptionsPattern[]] :=
Module[{tableIndex = 0, startTime = AbsoluteTime[]},
	Monitor[
		Table[tableIndex++;expr, n]
		,
		If[OptionValue[ProgressTable, {opts}, "ShowInfo"],
			DetailedIndicator[tableIndex, n, startTime, OptionValue[ProgressTable, {opts}, "Label"]]
			,
			DefaultIndicator[tableIndex, n]
		]
	]
];

ProgressTable[expr_, iterators__, opts:OptionsPattern[]]:= 
Module[{tuplesLength, tableIndex = 0, startTime = AbsoluteTime[]},
	(* This is beautiful. Author: Giovanni F. From: https://mathematica.stackexchange.com/a/86059 *)
	tuplesLength = Length[Tuples[Map[Range[Rest[#] /. List->Sequence]&, {iterators}]]];

	Monitor[
		Table[tableIndex++;expr, iterators]
		,
		If[OptionValue[ProgressTable, {opts}, "ShowInfo"],
			DetailedIndicator[tableIndex, tuplesLength, startTime, OptionValue[ProgressTable, {opts}, "Label"]]
			,
			DefaultIndicator[tableIndex, tuplesLength]
		]
	]
];


(* ::Section::Closed:: *)
(*Parallel tables with progress*)


SetAttributes[ProgressParallelTable,HoldFirst];
Options[ProgressParallelTable]={"ShowInfo"->True, "Label"->"Evaluating..."};

ProgressParallelTable[expr_, n_Integer, opts:OptionsPattern[{ProgressParallelTable, ParallelTable}]]:=
Module[{tableIndex = 0,startTime = AbsoluteTime[]},
	SetSharedVariable[tableIndex];

	Monitor[
		ParallelTable[tableIndex++;expr, n, Evaluate[FilterRules[{opts}, Options[ParallelTable]]]]
		,
		If[OptionValue[ProgressParallelTable, {opts}, "ShowInfo"],
			DetailedIndicator[tableIndex, n, startTime, OptionValue[ProgressParallelTable, {opts}, "Label"]]
			,
			DefaultIndicator[tableIndex, n]
		]
	]
];

ProgressParallelTable[expr_, iterators__, opts:OptionsPattern[{ProgressParallelTable, ParallelTable}]]:= 
Module[{tuplesLength, tableIndex = 0, startTime = AbsoluteTime[]},
	SetSharedVariable[tableIndex];
	tuplesLength = Length[Tuples[Map[Range[Rest[#] /. List->Sequence]&, {iterators}]]];

	Monitor[
		ParallelTable[tableIndex++;expr, iterators, Evaluate[FilterRules[{opts}, Options[ParallelTable]]]]
		,
		If[OptionValue[ProgressParallelTable, {opts}, "ShowInfo"],
			DetailedIndicator[tableIndex, tuplesLength, startTime, OptionValue[ProgressParallelTable, {opts}, "Label"]]
			,
			DefaultIndicator[tableIndex, tuplesLength]
		]
	]
];


(* ::Section::Closed:: *)
(*Nest with progress*)


ProgressNest[f_, expr_, n_, opts:OptionsPattern[]]:= Module[{nestIndex = 0, startTime = AbsoluteTime[]},
	Monitor[
		Nest[(nestIndex++; f[#])&, expr, n]
		,
		If[OptionValue[ProgressTable, {opts}, "ShowInfo"],
			DetailedIndicator[nestIndex, n, startTime, OptionValue[ProgressTable,{opts},"Label"]],
			DefaultIndicator[nestIndex, n]
		]
	]
];
ProgressNestList[f_, expr_, n_, opts:OptionsPattern[]] := Module[{nestIndex = 0, startTime = AbsoluteTime[]},
	Monitor[
		NestList[(nestIndex++;f[#])&,expr,n]
		,
		If[OptionValue[ProgressTable,{opts},"ShowInfo"],
			DetailedIndicator[nestIndex,n,startTime,OptionValue[ProgressTable,{opts},"Label"]],
			DefaultIndicator[nestIndex,n]
		]
	]
];


(* ::Chapter:: *)
(*End of package*)


End[ ]

EndPackage[ ]
