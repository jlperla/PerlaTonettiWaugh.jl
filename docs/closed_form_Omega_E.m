(* ::Package:: *)

(* This notebook implements the algebra to ensure that the boundary \
conditions are met for the ODE, and that E(t) = \[Delta]*)

ConvertToPureFunction[expr_, vars_List] :=
 With[{variables = Unevaluated@vars},
  Block[variables,
   Evaluate@(Hold[expr] /.
        Thread[vars -> Slot /@ Range@Length@vars]) & // ReleaseHold]]


(* LINEAR *)
$Assumptions = \[CapitalOmega]0 > 0 && \[CapitalOmega]T >
0 && \[CapitalOmega]T < \[CapitalOmega]0 && \[Delta] > 0 &&
T > 0 && a > 0 && b > 0;
Esub = {e -> (b # + a &)};(* Not scaled *)
\[CapitalOmega]sol = DSolve[{D[\[CapitalOmega][t], t] == (e[t] - \[Delta]) \[CapitalOmega][t], \[CapitalOmega][0] == \[CapitalOmega]0} /. Esub, \[CapitalOmega][t], t] //FullSimplify;
\[CapitalOmega]sub = {\[CapitalOmega] ->
ConvertToPureFunction[\[CapitalOmega]sol[[1, 1, 2]], {t}]};
\[CapitalOmega]terminal =
Log[\[CapitalOmega]T] == Log[\[CapitalOmega][T]] /.
Union[\[CapitalOmega]sub, Esub] // PowerExpand ;
coefficientssol = First@Solve[{\[Delta] == e[T] ,
   \[CapitalOmega]terminal} /.
  Union[\[CapitalOmega]sub, Esub], {a, b}] // FullSimplify;
Print["Linear Solution:"]
solutions =
Union[First@\[CapitalOmega]sol /. coefficientssol // PowerExpand //
FullSimplify,
{e[t] -> (e[t] //. Union[Esub, coefficientssol] // PowerExpand //
  FullSimplify)}]
puresolutions = {e -> (ConvertToPureFunction[
  solutions[[1,
    2]], {t}]), \[CapitalOmega] -> (ConvertToPureFunction[
  solutions[[2, 2]], {t}])};
  
Print["A Few Checks:"]
e[T] == \[Delta] /. puresolutions
\[CapitalOmega]0 == \[CapitalOmega][0] /. puresolutions
\[CapitalOmega]T == \[CapitalOmega][T] /. puresolutions

Print["Monotonicity?"]
e'[t]  > 0/. puresolutions //FullSimplify


(* QUADRATIC *)
$Assumptions = \[CapitalOmega]0 > 0 && \[CapitalOmega]T >
0 && \[CapitalOmega]T < \[CapitalOmega]0 && \[Delta] > 0 &&
T > 0 && a !=  0 && b !=  0 && c1 != 0 && Log[\[CapitalOmega]T]<Log[\[CapitalOmega]0];
Esub = {e -> (c1 #^2 + b # + a &)};(* Not scaled *)
\[CapitalOmega]sol =
DSolve[{D[\[CapitalOmega][t],
   t] == (e[t] - \[Delta]) \[CapitalOmega][t], \[CapitalOmega][
   0] == \[CapitalOmega]0} /. Esub, \[CapitalOmega][t], t] //
FullSimplify;
\[CapitalOmega]sub = {\[CapitalOmega] ->
ConvertToPureFunction[\[CapitalOmega]sol[[1, 1, 2]], {t}]};
\[CapitalOmega]terminal =
Log[\[CapitalOmega]T] == Log[\[CapitalOmega][T]] /.
Union[\[CapitalOmega]sub, Esub] // PowerExpand ;
coefficientssol = First@Solve[{\[Delta] == e[T] ,
  \[CapitalOmega]terminal} /. Union[\[CapitalOmega]sub, Esub], {a,
  b}] // FullSimplify
Print["Quadratic Solution:"]
solutions =
Union[First@\[CapitalOmega]sol /. coefficientssol // PowerExpand //
FullSimplify,
{e[t] -> (e[t] //. Union[Esub, coefficientssol] // PowerExpand //
  FullSimplify)}]
puresolutions = {e -> (ConvertToPureFunction[
  solutions[[1,
    2]], {t}]), \[CapitalOmega] -> (ConvertToPureFunction[
  solutions[[2, 2]], {t}])};
Print["A Few Checks:"]
e[T] == \[Delta] /. puresolutions
\[CapitalOmega]0 == \[CapitalOmega][0] /. puresolutions
\[CapitalOmega]T == \[CapitalOmega][T] /. puresolutions

Print["Monotonicity Check"]
mononticitycheck = e'[t]  > 0/. puresolutions //FullSimplify;
Reduce[(mononticitycheck /. t-> 0)&&(mononticitycheck /. t-> T)] //FullSimplify


(* CUBIC *)
$Assumptions = \[CapitalOmega]0 > 0 && \[CapitalOmega]T >
0 && \[CapitalOmega]T < \[CapitalOmega]0 && \[Delta] > 0 &&
T > 0 && a !=  0 && b !=  0 && c1 != 0 && c2 != 0 && Log[\[CapitalOmega]T]<Log[\[CapitalOmega]0];
Esub = {e -> (c2 #^3 + c1 #^2 + b # + a &)};(* Not scaled *)
\[CapitalOmega]sol =
DSolve[{D[\[CapitalOmega][t],
   t] == (e[t] - \[Delta]) \[CapitalOmega][t], \[CapitalOmega][
   0] == \[CapitalOmega]0} /. Esub, \[CapitalOmega][t], t] //
FullSimplify;
\[CapitalOmega]sub = {\[CapitalOmega] ->
ConvertToPureFunction[\[CapitalOmega]sol[[1, 1, 2]], {t}]};
\[CapitalOmega]terminal =
Log[\[CapitalOmega]T] == Log[\[CapitalOmega][T]] /.
Union[\[CapitalOmega]sub, Esub] // PowerExpand ;
coefficientssol = First@Solve[{\[Delta] == e[T] ,
  \[CapitalOmega]terminal} /. Union[\[CapitalOmega]sub, Esub], {a,
  b}] // FullSimplify
Print["Cubic Solution:"]
solutions =
Union[First@\[CapitalOmega]sol /. coefficientssol // PowerExpand //
FullSimplify,
{e[t] -> (e[t] //. Union[Esub, coefficientssol] // PowerExpand //
  FullSimplify)}]
puresolutions = {e -> (ConvertToPureFunction[
  solutions[[1,
    2]], {t}]), \[CapitalOmega] -> (ConvertToPureFunction[
  solutions[[2, 2]], {t}])};
Print["A Few Checks:"]
e[T] == \[Delta] /. puresolutions //FullSimplify
\[CapitalOmega]0 == \[CapitalOmega][0] /. puresolutions//FullSimplify
\[CapitalOmega]T == \[CapitalOmega][T] /. puresolutions//FullSimplify

Print["Monotonicity?"]
mononticitycheck = e[t] > 0 /. puresolutions //FullSimplify;
Reduce[(mononticitycheck /. t-> 0)&&(mononticitycheck /. t-> T)&&(mononticitycheck /. t-> T/2)] //FullSimplify


(* QUARTIC *)
$Assumptions = \[CapitalOmega]0 > 0 && \[CapitalOmega]T >
0 && \[CapitalOmega]T < \[CapitalOmega]0 && \[Delta] > 0 &&
T > 0 && a > 0 && b > 0 && c1 != 0 && c2 \[Element]Reals && c3 \[Element]Reals && Log[\[CapitalOmega]T]<Log[\[CapitalOmega]0];
Esub = {e -> (c3 #^4 + c2 #^3 + c1 #^2 + b # +
  a &)}; \[CapitalOmega]sol =
DSolve[{D[\[CapitalOmega][t],
   t] == (e[t] - \[Delta]) \[CapitalOmega][t], \[CapitalOmega][
   0] == \[CapitalOmega]0} /. Esub, \[CapitalOmega][t], t] //
FullSimplify;
\[CapitalOmega]sub = {\[CapitalOmega] ->
ConvertToPureFunction[\[CapitalOmega]sol[[1, 1, 2]], {t}]};
\[CapitalOmega]terminal =
Log[\[CapitalOmega]T] == Log[\[CapitalOmega][T]] /.
Union[\[CapitalOmega]sub, Esub] // PowerExpand ;
coefficientssol = First@Solve[{\[Delta] == e[T] ,
   \[CapitalOmega]terminal} /.
  Union[\[CapitalOmega]sub, Esub], {a, b}] // FullSimplify;
Print["Quartic Solution:"]
solutions =
Union[First@\[CapitalOmega]sol /. coefficientssol // PowerExpand //
FullSimplify,
{e[t] -> (e[t] //. Union[Esub, coefficientssol] // PowerExpand //
  FullSimplify)}]
  puresolutions = {e -> (ConvertToPureFunction[
  solutions[[1,
    2]], {t}]), \[CapitalOmega] -> (ConvertToPureFunction[
  solutions[[2, 2]], {t}])};
Print["A Few Checks:"]
e[T] == \[Delta] /. puresolutions //FullSimplify
\[CapitalOmega]0 == \[CapitalOmega][0] /. puresolutions//FullSimplify
\[CapitalOmega]T == \[CapitalOmega][T] /. puresolutions//FullSimplify

Print["Monotonicity?"]
mononticitycheck = e[t] > 0 /. puresolutions //FullSimplify;
Reduce[(mononticitycheck /. t-> 0)&&(mononticitycheck /. t-> T)&&(mononticitycheck /. t-> T/2)] //FullSimplify


(* QUINTIC *)
$Assumptions = \[CapitalOmega]0 >0 && \[CapitalOmega]T > 0 && \[CapitalOmega]T < \[CapitalOmega]0 && \[Delta] > 0 && T > 0 && a > 0 && b > 0 && c1 != 0 && c2 != 0 && c3 != 0 && c4 != 0 && Log[\[CapitalOmega]T]<Log[\[CapitalOmega]0];
Esub = {e -> (c4 #^5 + c3 #^4 + c2 #^3 + c1 #^2 + b # + a&)};\[CapitalOmega]sol = DSolve[{D[\[CapitalOmega][t], t]== (e[t] - \[Delta])\[CapitalOmega][t], \[CapitalOmega][0] == \[CapitalOmega]0} /. Esub, \[CapitalOmega][t],t] //FullSimplify;
\[CapitalOmega]sub = {\[CapitalOmega] -> ConvertToPureFunction[\[CapitalOmega]sol[[1,1,2]], {t}]};
\[CapitalOmega]terminal = Log[\[CapitalOmega]T]== Log[\[CapitalOmega][T]]/. Union[\[CapitalOmega]sub,Esub] //PowerExpand ;
coefficientssol = First@Solve[{\[Delta] == e[T] ,
\[CapitalOmega]terminal} /. Union[\[CapitalOmega]sub,Esub], {a,b}] //FullSimplify;
Print["Quintic Solution:"]
solutions = Union[First@\[CapitalOmega]sol /. coefficientssol //Simplify,
{e[t] -> (e[t] //. Union[Esub,coefficientssol])}//Simplify]
puresolutions = {e -> (ConvertToPureFunction[
  solutions[[1,
    2]], {t}]), \[CapitalOmega] -> (ConvertToPureFunction[
  solutions[[2, 2]], {t}])};
Print["A Few Checks:"]
e[T] == \[Delta] /. puresolutions //FullSimplify
\[CapitalOmega]0 == \[CapitalOmega][0] /. puresolutions//FullSimplify
\[CapitalOmega]T == \[CapitalOmega][T] /. puresolutions//FullSimplify


(* SEXTIC *)
$Assumptions = \[CapitalOmega]0 >0 && \[CapitalOmega]T > 0 && \[CapitalOmega]T < \[CapitalOmega]0 && \[Delta] > 0 && T > 0 && a > 0 && b > 0 && c1 != 0 && c2 != 0 && c3 != 0 && c4 != 0;
Esub = {e -> (c5 #^6 + c4 #^5 + c3 #^4 + c2 #^3 + c1 #^2 + b # + a&)};\[CapitalOmega]sol = DSolve[{D[\[CapitalOmega][t], t]== (e[t] - \[Delta])\[CapitalOmega][t], \[CapitalOmega][0] == \[CapitalOmega]0} /. Esub, \[CapitalOmega][t],t] //FullSimplify;
\[CapitalOmega]sub = {\[CapitalOmega] -> ConvertToPureFunction[\[CapitalOmega]sol[[1,1,2]], {t}]};
\[CapitalOmega]terminal = Log[\[CapitalOmega]T]== Log[\[CapitalOmega][T]]/. Union[\[CapitalOmega]sub,Esub] //PowerExpand ;
coefficientssol = First@Solve[{\[Delta] == e[T] ,
\[CapitalOmega]terminal} /. Union[\[CapitalOmega]sub,Esub], {a,b}] //FullSimplify;
Print["Sextic Solution:"]
solutions = Union[First@\[CapitalOmega]sol /. coefficientssol //Simplify,
{e[t] -> (e[t] //. Union[Esub,coefficientssol])}//Simplify]
puresolutions = {e -> (ConvertToPureFunction[
  solutions[[1,
    2]], {t}]), \[CapitalOmega] -> (ConvertToPureFunction[
  solutions[[2, 2]], {t}])};
Print["A Few Checks:"]
e[T] == \[Delta] /.puresolutions
\[CapitalOmega]0 == \[CapitalOmega][0]/. puresolutions
\[CapitalOmega]T == \[CapitalOmega][T]/. puresolutions


(* SEPTIC *)
$Assumptions = \[CapitalOmega]0 >0 && \[CapitalOmega]T > 0 && \[CapitalOmega]T < \[CapitalOmega]0 && \[Delta] > 0 && T > 0 && a > 0 && b > 0 && c1 != 0 && c2 != 0 && c3 != 0 && c4 != 0;
Esub = {e -> (c6 #^7 + c5 #^6 + c4 #^5 + c3 #^4 + c2 #^3 + c1 #^2 + b # + a&)};\[CapitalOmega]sol = DSolve[{D[\[CapitalOmega][t], t]== (e[t] - \[Delta])\[CapitalOmega][t], \[CapitalOmega][0] == \[CapitalOmega]0} /. Esub, \[CapitalOmega][t],t] //FullSimplify;
\[CapitalOmega]sub = {\[CapitalOmega] -> ConvertToPureFunction[\[CapitalOmega]sol[[1,1,2]], {t}]};
\[CapitalOmega]terminal = Log[\[CapitalOmega]T]== Log[\[CapitalOmega][T]]/. Union[\[CapitalOmega]sub,Esub] //PowerExpand ;
coefficientssol = First@Solve[{\[Delta] == e[T] ,
\[CapitalOmega]terminal} /. Union[\[CapitalOmega]sub,Esub], {a,b}] //FullSimplify;
Print["Septic Solution:"]
solutions = Union[First@\[CapitalOmega]sol /. coefficientssol //Simplify,
{e[t] -> (e[t] //. Union[Esub,coefficientssol])}//Simplify]
Print["A Few Checks:"]
e[T] == \[Delta] /.puresolutions
\[CapitalOmega]0 == \[CapitalOmega][0]/. puresolutions
\[CapitalOmega]T == \[CapitalOmega][T]/. puresolutions



