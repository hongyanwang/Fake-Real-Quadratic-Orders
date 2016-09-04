# Fake-Real-Quadratic-Orders
Computations for the Cohen-Lenstra Heuristic and the Ankeny-Artin-Chowla conjecture on Fake Real Quadratic Orders.

The program requires Sayles's libraries liboptarith and libqform
https://github.com/maxwellsayles

The computation is based on Mosunov's data: 
class numbers of imaginary quadratic fields with fundamental discriminants up to 2^40
http://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups

main.c is for the combination of the Cohen-Lenstra Heuristic and the Ankeny-Artin-Chowla conjecture with small prime set.
aacD.c is to look for counterexamples of the Ankeny-Artin-Chowla conjecture for selected D with p in a range.

More descriptions can be found in FRQO.pdf
