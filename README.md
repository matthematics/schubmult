# schubmult

## Program and package for rapid computation of Littlewood-Richardson coefficients of Schubert 
polynomials, with optional Sage integration

The main purpose of this python package is for executing scripts to compute coefficients of products of 
various types of Schubert polynomials. Coproducts can also be computed, as well as substitution of 
commuting difference operators for quantum double Schubert polynomials. Quantum multiplication also has 
parabolic subgroup support, computed via the Peterson-Woodward comparison theorem. **Note that except 
for quantum Schubert polynomial multiplication with the --basic-pieri option, the methodology for 
quantum/quantum double Schubert polynomials is conjectural at this time.**


[Docs to be hosted on Wiki](https://github.com/matthematics/schubmult/wiki/schubmult-home)


## Basic script command lines, one-line notation

### schubmult_py - ordinary Schubert polynmoials

```
usage: schubmult_py [-h] [-np] [--code] [--mult MULT [MULT ...]] [--coprod] [--display-mode 
{basic,pretty,latex,raw}] 
hyphen-separated list of perms [hyphen-separated list of perms ...]


Compute products of ordinary Schubert polynomials

positional arguments:
  hyphen-separated list of perms
                        Space-delimited permutations separated by hyphens, e. g. 3 4 1 2 - 5 1 2 4 3


options:
  -h, --help            show this help message and exit
  -np, --no-print       Compute the result but do not print it

  --code                Permutations represented by the Lehmer code

  --mult MULT [MULT ...]
                        Some additional terms in the ring to multiply by

  --coprod              Compute the coproduct (different syntax: one permutation, then a hyphen, the 
variable index 
positions to split on)
  --display-mode {basic,raw}
                        Method of displaying the output. Default basic


Example:
        schubmult_py 5 1 7 3 2 6 4 - 2 1 6 3 5 4
            or equivalently
        schubmult_py --code 4 0 4 1 0 1 - 1 0 3 0 1
            or alternatively
        schubmult_py --coprod --code 2 0 3 0 1 - 2 4
```


<!-- ## Quantum commuting difference operators

schubmult_q_double has a feature for displaying the coefficients of the divided difference operators in 
the evaluation of the quantum double Schubert polynomials on the commuting difference operators of 
Fomin, Gelfand, and Postnikov. It is necessary to cap the value of n in the group S_n we are working in 
because as n increases the expression does not stabilize.
```bash
schubmult_q_double --nil-hecke 6 --code 2 2 --display-positive

``` -->

### schubmult_double - Double Schubert polynmoials

```
usage: schubmult_double [-h] [-np] [--code] [--mult MULT [MULT ...]] [--coprod] [--display-positive] [--
optimizer-message] [--down] [-nc] 
[--mixed-var] [--expand] [--display-mode {basic,pretty,latex,raw}]

hyphen-separated list of perms [hyphen-separated list of perms ...]


Compute coefficients of product of double Schubert polynomials in the same or different sets of 
coefficient variables

positional arguments:
  hyphen-separated list of perms
                        Space-delimited permutations separated by hyphens, e. g. 3 4 1 2 - 5 1 2 4 3


options:
  -h, --help            show this help message and exit
  -np, --no-print       Compute the result but do not print it

  --code                Permutations represented by the Lehmer code

  --mult MULT [MULT ...]
                        Some additional terms in the ring to multiply by

  --coprod              Compute the coproduct (different syntax: one permutation, then a hyphen, the 
variable index 
  positions to split on)
  --display-positive    Display the result in terms of the positive roots, or if mixed variable attempt 
to display the 
  result as a positive algebraic combination of terms of the form y_i - z_j

  --optimizer-message   Display debug output during integer optimization for --display-positive

  --down                Reverse multiplication
  -nc, --no-check       Do not check if positive result matches the original

  --mixed-var           Used mixed variables y and z
  --expand              Expand the output rather than leaving it as originally computed (slow)

  --display-mode {basic,pretty,latex,raw}
                        Method of displaying the output. Default basic


Example:
        schubmult_double 5 1 7 3 2 6 4 - 2 1 6 3 5 4 [ --display-positive]

            or equivalently
        schubmult_double --code 4 0 4 1 0 1 - 1 0 3 0 1 [ --display-positive]

            or alternatively
        schubmult_double --coprod --code 2 0 3 0 1 - 2 4  [ --display-positive]

```

### schubmult_q - Quantum Schubert polynomials

```
usage: schubmult_q [-h] [-np] [--code] [--mult MULT [MULT ...]] [--parabolic PARABOLIC [PARABOLIC ...]] 
[--basic-pieri] 
[--display-mode {basic,pretty,latex,raw}]
                   hyphen-separated list of perms [hyphen-separated list of perms ...]


Compute products of quantum Schubert polynomials

positional arguments:
  hyphen-separated list of perms
                        Space-delimited permutations separated by hyphens, e. g. 3 4 1 2 - 5 1 2 4 3


options:
  -h, --help            show this help message and exit
  -np, --no-print       Compute the result but do not print it

  --code                Permutations represented by the Lehmer code

  --mult MULT [MULT ...]
                        Some additional terms in the ring to multiply by

  --parabolic PARABOLIC [PARABOLIC ...]
                        Generators of the parabolic subgroup to compute quantum coeffs for

  --basic-pieri         Do not apply conjectural computation optimization to quantum

  --display-mode {basic,pretty,latex,raw}
                        Method of displaying the output. Default basic


Example:
        schubmult_q 5 1 7 3 2 6 4 - 2 1 6 3 5 4
            or equivalently
        schubmult_q --code 4 0 4 1 0 1 - 1 0 3 0 1
```

### schubmult_q_double - Quantum double Schubert polynomials


```
usage: schubmult_q_double [-h] [-np] [--code] [--mult MULT [MULT ...]] [--display-positive] [--optimizer-
message] [--down] 
[-nc] [--mixed-var] [--expand] [--parabolic PARABOLIC [PARABOLIC ...]] [--basic-pieri]

                          [--nil-hecke N] [--nil-hecke-apply N] [--display-mode 
{basic,pretty,latex,raw}]
                          hyphen-separated list of perms [hyphen-separated list of perms ...]


Compute coefficients of products of quantum double Schubert polynomials in the same or different sets of 
coefficient variables

positional arguments:
  hyphen-separated list of perms
                        Space-delimited permutations separated by hyphens, e. g. 3 4 1 2 - 5 1 2 4 3


options:
  -h, --help            show this help message and exit
  -np, --no-print       Compute the result but do not print it

  --code                Permutations represented by the Lehmer code

  --mult MULT [MULT ...]
                        Some additional terms in the ring to multiply by

  --display-positive    Display the result in terms of the positive roots, or if mixed variable attempt 
to display the 
result as a positive algebraic combination of terms of the form y_i - z_j

  --optimizer-message   Display debug output during integer optimization for --display-positive

  --down                Reverse multiplication
  -nc, --no-check       Do not check if positive result matches the original

  --mixed-var           Used mixed variables y and z
  --expand              Expand the output rather than leaving it as originally computed (slow)

  --parabolic PARABOLIC [PARABOLIC ...]
                        Generators of the parabolic subgroup to compute quantum coeffs for

  --basic-pieri         Do not apply conjectural computation optimization to quantum

  --nil-hecke N         Substitute up to N of Fomin-Gelfand-Postnikov commuting difference operators

  --nil-hecke-apply N   Substitute commuting difference operators for perm1, then apply to Schub indexed 
by perm2
  --display-mode {basic,pretty,latex,raw}
                        Method of displaying the output. Default basic


Example:
        schubmult_q_double 5 1 7 3 2 6 4 - 2 1 6 3 5 4 [ --display-positive]

            or equivalently
        schubmult_q_double --code 4 0 4 1 0 1 - 1 0 3 0 1 [ --display-positive]

```

## Diplaying the result positively

The command line argument `--display-positive `is available in schubmult_double and schubmult_q_double, 
which displays the result positively (if possible, this is still only always possible conjecturally). It 
will fail and print out the offending case if it finds a counterexample. This is highly processor 
intensive.


Runtime will vary tremendously by case. The general problem is #P-hard. Though the result is always 
nonnegative (which at least is known for schubmult_py, schubmult_q, schubmult_double, and 
schubmult_q_double) and the problem is in GapP, it is not known to be in #P at this time.


schubmult_py is for multiplying ordinary Schubert polynomials. schubmult_double is for multiplying 
double Schubert polynomials in different sets of coefficient variables (`--mixed-var`) or in the same 
set of coefficient variables (by default). Similarly, schubmult_q is for multiplying quantum Schubert 
polynomials, schubmult_q_double is for multiplying quantum double Schubert polynomials (in different 
sets of coefficient variables with the `--mixed-var` option, or the same set), or in other words it 
computes the Gromov-Witten invariants, equivariant Gromov-Witten invariants, and (mixed?) equivariant 
Gromov-Witten invariants of the complete flag variety. All have the same command line syntax as 
schubmult, except when using the --code option. schubmult_double/schubmult_q_double display the result 
with nonnegative coefficients (and the q variables) with the `--display-positive` option (either in the 
negative simple roots r_i, or in y_i - z_j for `--mixed-var`).


schubmult_xx --coprod allows you to split (double) Schubert polynomials along certain indices (not 
available for quantum). It takes one permutation as an argument, followed by a dash -, then the set of 
indices you would like to split on.

# Sage integration (as of version 2.0.0)

[SageMath](https://www.sagemath.org/) is a computer algebra system that, while wonderful, is monstrously 
large and only works on posix-based operating systems (including WSL VMs, so it is still usable on 
Windows). This is why Sage support is provided optionally in schubmult. The syntax to install the Sage 
dependencies is

```
pip install schubmult[sage]
```

This will install the [sagemath-standard](https://pypi.org/project/sagemath-standard/) python package in 
addition to the other dependencies. **Again, this only works on Linux, MacOS, or WSL.** To use with a 
currently installed SageMath distribution, use sage's python interpreter to install the package (the 
`[sage]` piece is not required in that case).

## Basic sage example

```python
sage: from schubmult.sage_integration import FastSchubertPolynomialRing, 
FastDoubleSchubertPolynomialRing,
FastQuantumSchubertPolynomialRing, FastQuantumDoubleSchubertPolynomialRing

sage: SingleRing = FastSchubertPolynomialRing(ZZ, 100, "x")
sage: SingleRing([3,4,1,2])
Sx([3, 4, 1, 2])
sage: SingleRing([3,4,1,2]) * SingleRing([5,1,4,2,3])
Sx([7, 3, 4, 1, 2, 5, 6]) + Sx([7, 4, 2, 1, 3, 5, 6]) + Sx([7, 5, 1, 2, 3, 4, 6])

```

## Mixed variable (Molev-Sagan) type products

```python
sage: DoubleRing = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
sage: DoubleRing([3,4,1,2])*DoubleRing([1,4,2,3])
(y1*y2-y1*y4-y2*y4+y4^2)*DSx([3, 4, 1, 2], 'y') + (-y1-y2+y4+y5)*DSx([3, 5, 1, 2, 4], 'y') + DSx([3, 6, 
1, 2, 4, 5], 'y')
sage: DoubleRing([3,4,1,2]) * DoubleRing([1,4,2,3],"z")
(y3^2+y3*y4+y4^2-y3*z1-y4*z1-y3*z2-y4*z2+z1*z2-y3*z3-y4*z3+z1*z3+z2*z3)*DSx([3, 4, 1, 2], 'y') + (y3+y4+
y5-z1-z2-z3)*DSx([3, 5, 1, 2, 4], 'y') + DSx([3, 6, 1, 2, 4, 5], 'y')
```

## expand()

```python
sage: SingleRingQ = FastQuantumSchubertPolynomialRing(ZZ, 100, "x")
sage: SingleRingQ([2,3,1,4]).expand()
x1*x2 + q_1
```

## Coercion

Coercion was implemented as widely as possible.
```python
sage: DoubleRing([1,4,2,3],"z") * SingleRing([3,4,1,2])
z1^2*z4^2*DSx([1, 4, 2, 3], 'z') + (z1^2*z4+z1^2*z5)*DSx([1, 5, 2, 3, 4], 'z') + z1^2*DSx([1, 6, 2, 3, 
4, 5], 'z') + (z1*z4^2+z2*z4^2)*DSx([2, 4, 1, 3], 'z') + (z1*z4+z2*z4+z1*z5+z2*z5)*DSx([2, 5, 1, 3, 4], 
'z') + (z1+z2)*DSx([2, 6, 1, 3, 4, 5], 'z') + z4^2*DSx([3, 4, 1, 2], 'z') + (z4+z5)*DSx([3, 5, 1, 2, 4], 
'z') + DSx([3, 6, 1, 2, 4, 5], 'z')
sage: SingleRingQ([2,3,1,4]) * SingleRing([4,1,3,2])
(-2*q1^2*q2+q1*q2*q3)*QSx([1]) + q1*q2*QSx([1, 3, 4, 2]) + (-q1^2)*QSx([2, 3, 1]) + q1*QSx([2, 4, 3, 1]) 
+ (-q1*q2)*QSx([3, 1, 2]) + q1*QSx([3, 2, 4, 1]) + q1*QSx([3, 4, 1, 2]) + (-q1)*QSx([4, 2, 1, 3]) + 
QSx([5, 2, 3, 1, 4]) + QSx([5, 3, 1, 2, 4])
sage: R.<x1, x2> = PolynomialRing(ZZ, 2)
sage: SingleRing([1,3,2]) - x1 - x2 == 0
True
```

## Coproducts

FastSchubertPolynomialRing and FastDoubleSchubertPolynomialRings are bialgebras and each element 
implements the `coproduct()` member function. `set_coproduct_indices()` on the base ring will determine 
the variables to partition on.
```ada
sage: DoubleRing.set_coproduct_indices((1,3))
sage: DoubleRing([4,1,5,2,3], "z").coproduct()
(y1^2-y1*z2-y1*z3+z2*z3)*DSx([4, 1, 2, 3], 'z') # DSx([1], 'y') + (y1+y2-z2-z3)*DSx([4, 1, 2, 3], 'z') # 
DSx([2, 1], 'y') + DSx([4, 1, 2, 3], 'z') # DSx([3, 1, 2], 'y') + (y1-z3)*DSx([4, 2, 1, 3], 'z') # 
DSx([1], 'y') + DSx([4, 2, 1, 3], 'z') # DSx([2, 1], 'y') + DSx([4, 3, 1, 2], 'z') # DSx([1], 'y')

```

## Demonstration of quantum double mixed products

```python
sage: QuantumDoubleRing = FastQuantumDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y","z"))
sage: QuantumDoubleRing([4,1,3,2])*QuantumDoubleRing([5,1,3,2,4], "z")
(q1*q2*q3*y1^3+q1*q2*q3*y1^2*y4+q1*q2*q3*y1*y4^2+q1*q2*q3*y4^3-q1*q2*q3*y1^2*z1-q1*q2*q3*y1*y4*z1-
q1*q2*q3*y4^2*z1-q1*q2*q3*y1^2*z2-q1*q2*q3*y1*y4*z2-q1*q2*q3*y4^2*z2+q1*q2*q3*y1*z1*z2+q1*q2*q3*y4*z1*z2-
q1*q2*q3*y1^2*z3-q1*q2*q3*y1*y4*z3-q1*q2*q3*y4^2*z3+q1*q2*q3*y1*z1*z3+q1*q2*q3*y4*z1*z3+
q1*q2*q3*y1*z2*z3+q1*q2*q3*y4*z2*z3-q1*q2*q3*z1*z2*z3-q1*q2*q3*y1^2*z4-q1*q2*q3*y1*y4*z4-
q1*q2*q3*y4^2*z4+q1*q2*q3*y1*z1*z4+q1*q2*q3*y4*z1*z4+q1*q2*q3*y1*z2*z4+q1*q2*q3*y4*z2*z4-
q1*q2*q3*z1*z2*z4+q1*q2*q3*y1*z3*z4+q1*q2*q3*y4*z3*z4-q1*q2*q3*z1*z3*z4-q1*q2*q3*z2*z3*z4)*QDSx([1], 
'y') + (q1*q2*q3*y1^2+q1*q2*q3*y1*y4+q1*q2*q3*y4^2+q1*q2*q3*y1*y5+q1*q2*q3*y4*y5+q1*q2*q3*y5^2-
q1*q2*q3*y1*z1-q1*q2*q3*y4*z1-q1*q2*q3*y5*z1-q1*q2*q3*y1*z2-q1*q2*q3*y4*z2-q1*q2*q3*y5*z2+q1*q2*q3*z1*z2-
q1*q2*q3*y1*z3-q1*q2*q3*y4*z3-q1*q2*q3*y5*z3+q1*q2*q3*z1*z3+q1*q2*q3*z2*z3-q1*q2*q3*y1*z4-q1*q2*q3*y4*z4-
q1*q2*q3*y5*z4+q1*q2*q3*z1*z4+q1*q2*q3*z2*z4+q1*q2*q3*z3*z4)*QDSx([1, 2, 3, 5, 4], 'y') + (q1*q2*q3*y1+
q1*q2*q3*y4+q1*q2*q3*y5+q1*q2*q3*y6-q1*q2*q3*z1-q1*q2*q3*z2-q1*q2*q3*z3-q1*q2*q3*z4)*QDSx([1, 2, 3, 6, 
4, 5], 'y') + q1*q2*q3*QDSx([1, 2, 3, 7, 4, 5, 6], 'y') + (q1*q3*y1^3+q1*q3*y1^2*y4+q1*q3*y1*y4^2+
q1*q3*y4^3-q1*q3*y1^2*z1-q1*q3*y1*y4*z1-q1*q3*y4^2*z1-q1*q3*y1^2*z2-q1*q3*y1*y4*z2-q1*q3*y4^2*z2+
q1*q3*y1*z1*z2+q1*q3*y4*z1*z2-q1*q3*y1^2*z3-q1*q3*y1*y4*z3-q1*q3*y4^2*z3+q1*q3*y1*z1*z3+q1*q3*y4*z1*z3+
q1*q3*y1*z2*z3+q1*q3*y4*z2*z3-q1*q3*z1*z2*z3-q1*q3*y1^2*z4-q1*q3*y1*y4*z4-q1*q3*y4^2*z4+q1*q3*y1*z1*z4+
q1*q3*y4*z1*z4+q1*q3*y1*z2*z4+q1*q3*y4*z2*z4-q1*q3*z1*z2*z4+q1*q3*y1*z3*z4+q1*q3*y4*z3*z4-q1*q3*z1*z3*z4-
q1*q3*z2*z3*z4)*QDSx([1, 4, 2, 3], 'y') + (q1*y1^3*y3+q1*y1^3*y4+q1*y1^2*y3*y4+q1*y1^2*y4^2+
q1*y1*y3*y4^2+q1*y1*y4^3+q1*y3*y4^3-q1*y1^3*z1-q1*y1^2*y3*z1-2*q1*y1^2*y4*z1-q1*y1*y3*y4*z1-
2*q1*y1*y4^2*z1-q1*y3*y4^2*z1-q1*y4^3*z1+q1*y1^2*z1^2+q1*y1*y4*z1^2+q1*y4^2*z1^2-q1*y1^3*z2-
q1*y1^2*y3*z2-2*q1*y1^2*y4*z2-q1*y1*y3*y4*z2-2*q1*y1*y4^2*z2-q1*y3*y4^2*z2-q1*y4^3*z2+2*q1*y1^2*z1*z2+
q1*y1*y3*z1*z2+3*q1*y1*y4*z1*z2+q1*y3*y4*z1*z2+2*q1*y4^2*z1*z2-q1*y1*z1^2*z2-q1*y4*z1^2*z2+q1*y1^2*z2^2+
q1*y1*y4*z2^2+q1*y4^2*z2^2-q1*y1*z1*z2^2-q1*y4*z1*z2^2-q1*y1^2*y3*z3-q1*y1^2*y4*z3-q1*y1*y3*y4*z3-
q1*y1*y4^2*z3-q1*y3*y4^2*z3+q1*y1^2*z1*z3+q1*y1*y3*z1*z3+2*q1*y1*y4*z1*z3+q1*y3*y4*z1*z3+q1*y4^2*z1*z3-
q1*y1*z1^2*z3-q1*y4*z1^2*z3+q1*y1^2*z2*z3+q1*y1*y3*z2*z3+2*q1*y1*y4*z2*z3+q1*y3*y4*z2*z3+q1*y4^2*z2*z3-
2*q1*y1*z1*z2*z3-q1*y3*z1*z2*z3-2*q1*y4*z1*z2*z3+q1*z1^2*z2*z3-q1*y1*z2^2*z3-q1*y4*z2^2*z3+q1*z1*z2^2*z3-
q1*y1^2*y3*z4-q1*y1^2*y4*z4-q1*y1*y3*y4*z4-q1*y1*y4^2*z4-q1*y3*y4^2*z4+q1*y1^2*z1*z4+q1*y1*y3*z1*z4+
2*q1*y1*y4*z1*z4+q1*y3*y4*z1*z4+q1*y4^2*z1*z4-q1*y1*z1^2*z4-q1*y4*z1^2*z4+q1*y1^2*z2*z4+q1*y1*y3*z2*z4+
2*q1*y1*y4*z2*z4+q1*y3*y4*z2*z4+q1*y4^2*z2*z4-2*q1*y1*z1*z2*z4-q1*y3*z1*z2*z4-2*q1*y4*z1*z2*z4+
q1*z1^2*z2*z4-q1*y1*z2^2*z4-q1*y4*z2^2*z4+q1*z1*z2^2*z4+q1*y1*y3*z3*z4+q1*y1*y4*z3*z4+q1*y3*y4*z3*z4-
q1*y1*z1*z3*z4-q1*y3*z1*z3*z4-q1*y4*z1*z3*z4+q1*z1^2*z3*z4-q1*y1*z2*z3*z4-q1*y3*z2*z3*z4-q1*y4*z2*z3*z4+
q1*z1*z2*z3*z4+q1*z2^2*z3*z4)*QDSx([1, 4, 3, 2], 'y') + (q1*y1^3+q1*y1^2*y4+q1*y1*y4^2+q1*y4^3-
q1*y1^2*z1-q1*y1*y4*z1-q1*y4^2*z1-q1*y1^2*z2-q1*y1*y4*z2-q1*y4^2*z2+q1*y1*z1*z2+q1*y4*z1*z2-q1*y1^2*z3-
q1*y1*y4*z3-q1*y4^2*z3+q1*y1*z1*z3+q1*y4*z1*z3+q1*y1*z2*z3+q1*y4*z2*z3-q1*z1*z2*z3-q1*y1^2*z4-
q1*y1*y4*z4-q1*y4^2*z4+q1*y1*z1*z4+q1*y4*z1*z4+q1*y1*z2*z4+q1*y4*z2*z4-q1*z1*z2*z4+q1*y1*z3*z4+
q1*y4*z3*z4-q1*z1*z3*z4-q1*z2*z3*z4)*QDSx([1, 4, 5, 2, 3], 'y') + (q1*q3*y1^2+q1*q3*y1*y4+q1*q3*y4^2+
q1*q3*y1*y5+q1*q3*y4*y5+q1*q3*y5^2-q1*q3*y1*z1-q1*q3*y4*z1-q1*q3*y5*z1-q1*q3*y1*z2-q1*q3*y4*z2-
q1*q3*y5*z2+q1*q3*z1*z2-q1*q3*y1*z3-q1*q3*y4*z3-q1*q3*y5*z3+q1*q3*z1*z3+q1*q3*z2*z3-q1*q3*y1*z4-
q1*q3*y4*z4-q1*q3*y5*z4+q1*q3*z1*z4+q1*q3*z2*z4+q1*q3*z3*z4)*QDSx([1, 5, 2, 3, 4], 'y') + (q1*y1^3+
q1*y1^2*y3+q1*y1^2*y4+q1*y1*y3*y4+q1*y1*y4^2+q1*y3*y4^2+q1*y1^2*y5+q1*y1*y3*y5+q1*y1*y4*y5+q1*y3*y4*y5+
q1*y1*y5^2+q1*y3*y5^2-2*q1*y1^2*z1-q1*y1*y3*z1-2*q1*y1*y4*z1-q1*y3*y4*z1-q1*y4^2*z1-2*q1*y1*y5*z1-
q1*y3*y5*z1-q1*y4*y5*z1-q1*y5^2*z1+q1*y1*z1^2+q1*y4*z1^2+q1*y5*z1^2-2*q1*y1^2*z2-q1*y1*y3*z2-
2*q1*y1*y4*z2-q1*y3*y4*z2-q1*y4^2*z2-2*q1*y1*y5*z2-q1*y3*y5*z2-q1*y4*y5*z2-q1*y5^2*z2+3*q1*y1*z1*z2+
q1*y3*z1*z2+2*q1*y4*z1*z2+2*q1*y5*z1*z2-q1*z1^2*z2+q1*y1*z2^2+q1*y4*z2^2+q1*y5*z2^2-q1*z1*z2^2-
q1*y1^2*z3-q1*y1*y3*z3-q1*y1*y4*z3-q1*y3*y4*z3-q1*y1*y5*z3-q1*y3*y5*z3+2*q1*y1*z1*z3+q1*y3*z1*z3+
q1*y4*z1*z3+q1*y5*z1*z3-q1*z1^2*z3+2*q1*y1*z2*z3+q1*y3*z2*z3+q1*y4*z2*z3+q1*y5*z2*z3-2*q1*z1*z2*z3-
q1*z2^2*z3-q1*y1^2*z4-q1*y1*y3*z4-q1*y1*y4*z4-q1*y3*y4*z4-q1*y1*y5*z4-q1*y3*y5*z4+2*q1*y1*z1*z4+
q1*y3*z1*z4+q1*y4*z1*z4+q1*y5*z1*z4-q1*z1^2*z4+2*q1*y1*z2*z4+q1*y3*z2*z4+q1*y4*z2*z4+q1*y5*z2*z4-
2*q1*z1*z2*z4-q1*z2^2*z4+q1*y1*z3*z4+q1*y3*z3*z4-q1*z1*z3*z4-q1*z2*z3*z4)*QDSx([1, 5, 3, 2, 4], 'y') + 
(q1*y1^2+q1*y1*y4+q1*y4^2+q1*y1*y5+q1*y4*y5+q1*y5^2-q1*y1*z1-q1*y4*z1-q1*y5*z1-q1*y1*z2-q1*y4*z2-
q1*y5*z2+q1*z1*z2-q1*y1*z3-q1*y4*z3-q1*y5*z3+q1*z1*z3+q1*z2*z3-q1*y1*z4-q1*y4*z4-q1*y5*z4+q1*z1*z4+
q1*z2*z4+q1*z3*z4)*QDSx([1, 5, 4, 2, 3], 'y') + (q1*q3*y1+q1*q3*y4+q1*q3*y5+q1*q3*y6-q1*q3*z1-q1*q3*z2-
q1*q3*z3-q1*q3*z4)*QDSx([1, 6, 2, 3, 4, 5], 'y') + (q1*y1^2+q1*y1*y3+q1*y1*y4+q1*y3*y4+q1*y1*y5+q1*y3*y5+
q1*y1*y6+q1*y3*y6-2*q1*y1*z1-q1*y3*z1-q1*y4*z1-q1*y5*z1-q1*y6*z1+q1*z1^2-2*q1*y1*z2-q1*y3*z2-q1*y4*z2-
q1*y5*z2-q1*y6*z2+2*q1*z1*z2+q1*z2^2-q1*y1*z3-q1*y3*z3+q1*z1*z3+q1*z2*z3-q1*y1*z4-q1*y3*z4+q1*z1*z4+
q1*z2*z4)*QDSx([1, 6, 3, 2, 4, 5], 'y') + (q1*y1+q1*y4+q1*y5+q1*y6-q1*z1-q1*z2-q1*z3-q1*z4)*QDSx([1, 6, 
4, 2, 3, 5], 'y') + q1*q3*QDSx([1, 7, 2, 3, 4, 5, 6], 'y') + (q1*y1+q1*y3-q1*z1-q1*z2)*QDSx([1, 7, 3, 2, 
4, 5, 6], 'y') + q1*QDSx([1, 7, 4, 2, 3, 5, 6], 'y') + (q1*q2*q3*y1^2+q1*q2*q3*y1*y2+q1*q2*q3*y2^2+
q1*q2*q3*y1*y4+q1*q2*q3*y2*y4+q1*q2*q3*y4^2-q1*q2*q3*y1*z1-q1*q2*q3*y2*z1-q1*q2*q3*y4*z1-q1*q2*q3*y1*z2-
q1*q2*q3*y2*z2-q1*q2*q3*y4*z2+q1*q2*q3*z1*z2-q1*q2*q3*y1*z3-q1*q2*q3*y2*z3-q1*q2*q3*y4*z3+q1*q2*q3*z1*z3+
q1*q2*q3*z2*z3-q1*q2*q3*y1*z4-q1*q2*q3*y2*z4-q1*q2*q3*y4*z4+q1*q2*q3*z1*z4+q1*q2*q3*z2*z4+
q1*q2*q3*z3*z4)*QDSx([2, 1], 'y') + (q1*q2*q3*y1+q1*q2*q3*y2+q1*q2*q3*y4+q1*q2*q3*y5-q1*q2*q3*z1-
q1*q2*q3*z2-q1*q2*q3*z3-q1*q2*q3*z4)*QDSx([2, 1, 3, 5, 4], 'y') + q1*q2*q3*QDSx([2, 1, 3, 6, 4, 5], 'y') 
+ (q1*q3*y1^2+q1*q3*y1*y2+q1*q3*y2^2+q1*q3*y1*y4+q1*q3*y2*y4+q1*q3*y4^2-q1*q3*y1*z1-q1*q3*y2*z1-
q1*q3*y4*z1-q1*q3*y1*z2-q1*q3*y2*z2-q1*q3*y4*z2+q1*q3*z1*z2-q1*q3*y1*z3-q1*q3*y2*z3-q1*q3*y4*z3+
q1*q3*z1*z3+q1*q3*z2*z3-q1*q3*y1*z4-q1*q3*y2*z4-q1*q3*y4*z4+q1*q3*z1*z4+q1*q3*z2*z4+
q1*q3*z3*z4)*QDSx([2, 4, 1, 3], 'y') + (q1*y1^2*y3+q1*y1*y2*y3+q1*y2^2*y3+q1*y1^2*y4+q1*y1*y2*y4+
q1*y2^2*y4+q1*y1*y3*y4+q1*y2*y3*y4+q1*y1*y4^2+q1*y2*y4^2+q1*y3*y4^2+q1*y4^3-q1*y1^2*z1-q1*y1*y2*z1-
q1*y2^2*z1-q1*y1*y3*z1-q1*y2*y3*z1-2*q1*y1*y4*z1-2*q1*y2*y4*z1-q1*y3*y4*z1-2*q1*y4^2*z1+q1*y1*z1^2+
q1*y2*z1^2+q1*y4*z1^2-q1*y1^2*z2-q1*y1*y2*z2-q1*y2^2*z2-q1*y1*y3*z2-q1*y2*y3*z2-2*q1*y1*y4*z2-
2*q1*y2*y4*z2-q1*y3*y4*z2-2*q1*y4^2*z2+2*q1*y1*z1*z2+2*q1*y2*z1*z2+q1*y3*z1*z2+3*q1*y4*z1*z2-q1*z1^2*z2+
q1*y1*z2^2+q1*y2*z2^2+q1*y4*z2^2-q1*z1*z2^2-q1*y1*y3*z3-q1*y2*y3*z3-q1*y1*y4*z3-q1*y2*y4*z3-q1*y3*y4*z3-
q1*y4^2*z3+q1*y1*z1*z3+q1*y2*z1*z3+q1*y3*z1*z3+2*q1*y4*z1*z3-q1*z1^2*z3+q1*y1*z2*z3+q1*y2*z2*z3+
q1*y3*z2*z3+2*q1*y4*z2*z3-2*q1*z1*z2*z3-q1*z2^2*z3-q1*y1*y3*z4-q1*y2*y3*z4-q1*y1*y4*z4-q1*y2*y4*z4-
q1*y3*y4*z4-q1*y4^2*z4+q1*y1*z1*z4+q1*y2*z1*z4+q1*y3*z1*z4+2*q1*y4*z1*z4-q1*z1^2*z4+q1*y1*z2*z4+
q1*y2*z2*z4+q1*y3*z2*z4+2*q1*y4*z2*z4-2*q1*z1*z2*z4-q1*z2^2*z4+q1*y3*z3*z4+q1*y4*z3*z4-q1*z1*z3*z4-
q1*z2*z3*z4)*QDSx([2, 4, 3, 1], 'y') + (q1*y1^2+q1*y1*y2+q1*y2^2+q1*y1*y4+q1*y2*y4+q1*y4^2-q1*y1*z1-
q1*y2*z1-q1*y4*z1-q1*y1*z2-q1*y2*z2-q1*y4*z2+q1*z1*z2-q1*y1*z3-q1*y2*z3-q1*y4*z3+q1*z1*z3+q1*z2*z3-
q1*y1*z4-q1*y2*z4-q1*y4*z4+q1*z1*z4+q1*z2*z4+q1*z3*z4)*QDSx([2, 4, 5, 1, 3], 'y') + (q1*q3*y1+q1*q3*y2+
q1*q3*y4+q1*q3*y5-q1*q3*z1-q1*q3*z2-q1*q3*z3-q1*q3*z4)*QDSx([2, 5, 1, 3, 4], 'y') + (q1*y1^2+q1*y1*y2+
q1*y2^2+q1*y1*y3+q1*y2*y3+q1*y1*y4+q1*y2*y4+q1*y3*y4+q1*y4^2+q1*y1*y5+q1*y2*y5+q1*y3*y5+q1*y4*y5+q1*y5^2-
2*q1*y1*z1-2*q1*y2*z1-q1*y3*z1-2*q1*y4*z1-2*q1*y5*z1+q1*z1^2-2*q1*y1*z2-2*q1*y2*z2-q1*y3*z2-2*q1*y4*z2-
2*q1*y5*z2+3*q1*z1*z2+q1*z2^2-q1*y1*z3-q1*y2*z3-q1*y3*z3-q1*y4*z3-q1*y5*z3+2*q1*z1*z3+2*q1*z2*z3-
q1*y1*z4-q1*y2*z4-q1*y3*z4-q1*y4*z4-q1*y5*z4+2*q1*z1*z4+2*q1*z2*z4+q1*z3*z4)*QDSx([2, 5, 3, 1, 4], 'y') +
 (q1*y1+q1*y2+q1*y4+q1*y5-q1*z1-q1*z2-q1*z3-q1*z4)*QDSx([2, 5, 4, 1, 3], 'y') + q1*q3*QDSx([2, 6, 1, 3, 
4, 5], 'y') + (q1*y1+q1*y2+q1*y3+q1*y4+q1*y5+q1*y6-2*q1*z1-2*q1*z2-q1*z3-q1*z4)*QDSx([2, 6, 3, 1, 4, 5], 
'y') + q1*QDSx([2, 6, 4, 1, 3, 5], 'y') + q1*QDSx([2, 7, 3, 1, 4, 5, 6], 'y') + (q1*q2*q3*y1+q1*q2*q3*y2+
q1*q2*q3*y3+q1*q2*q3*y4-q1*q2*q3*z1-q1*q2*q3*z2-q1*q2*q3*z3-q1*q2*q3*z4)*QDSx([3, 1, 2], 'y') + 
q1*q2*q3*QDSx([3, 1, 2, 5, 4], 'y') + (q1*y1^2*y3+q1*y1*y3^2+q1*y1^2*y4+2*q1*y1*y3*y4+q1*y3^2*y4+
q1*y1*y4^2+q1*y3*y4^2-q1*y1^2*z1-2*q1*y1*y3*z1-q1*y3^2*z1-2*q1*y1*y4*z1-2*q1*y3*y4*z1-q1*y4^2*z1+
q1*y1*z1^2+q1*y3*z1^2+q1*y4*z1^2-q1*y1^2*z2-2*q1*y1*y3*z2-q1*y3^2*z2-2*q1*y1*y4*z2-2*q1*y3*y4*z2-
q1*y4^2*z2+2*q1*y1*z1*z2+2*q1*y3*z1*z2+2*q1*y4*z1*z2-q1*z1^2*z2+q1*y1*z2^2+q1*y3*z2^2+q1*y4*z2^2-
q1*z1*z2^2-q1*y1*y3*z3-q1*y1*y4*z3-q1*y3*y4*z3+q1*y1*z1*z3+q1*y3*z1*z3+q1*y4*z1*z3-q1*z1^2*z3+
q1*y1*z2*z3+q1*y3*z2*z3+q1*y4*z2*z3-q1*z1*z2*z3-q1*z2^2*z3-q1*y1*y3*z4-q1*y1*y4*z4-q1*y3*y4*z4+
q1*y1*z1*z4+q1*y3*z1*z4+q1*y4*z1*z4-q1*z1^2*z4+q1*y1*z2*z4+q1*y3*z2*z4+q1*y4*z2*z4-q1*z1*z2*z4-
q1*z2^2*z4+q1*q3*y1+q1*q3*y2+q1*q3*y3+q1*q3*y4-q1*q3*z1-q1*q3*z2-q1*q3*z3-q1*q3*z4)*QDSx([3, 4, 1, 2], 
'y') + (q1*y1*y3+q1*y2*y3+q1*y3^2+q1*y1*y4+q1*y2*y4+2*q1*y3*y4+q1*y4^2-q1*y1*z1-q1*y2*z1-2*q1*y3*z1-
2*q1*y4*z1+q1*z1^2-q1*y1*z2-q1*y2*z2-2*q1*y3*z2-2*q1*y4*z2+2*q1*z1*z2+q1*z2^2-q1*y3*z3-q1*y4*z3+q1*z1*z3+
q1*z2*z3-q1*y3*z4-q1*y4*z4+q1*z1*z4+q1*z2*z4)*QDSx([3, 4, 2, 1], 'y') + (q1*y1+q1*y2+q1*y3+q1*y4-q1*z1-
q1*z2-q1*z3-q1*z4)*QDSx([3, 4, 5, 1, 2], 'y') + (q1*y1^2+2*q1*y1*y3+q1*y3^2+q1*y1*y4+q1*y3*y4+q1*y1*y5+
q1*y3*y5-2*q1*y1*z1-2*q1*y3*z1-q1*y4*z1-q1*y5*z1+q1*z1^2-2*q1*y1*z2-2*q1*y3*z2-q1*y4*z2-q1*y5*z2+
2*q1*z1*z2+q1*z2^2-q1*y1*z3-q1*y3*z3+q1*z1*z3+q1*z2*z3-q1*y1*z4-q1*y3*z4+q1*z1*z4+q1*z2*z4+
q1*q3)*QDSx([3, 5, 1, 2, 4], 'y') + (q1*y1+q1*y2+2*q1*y3+q1*y4+q1*y5-2*q1*z1-2*q1*z2-q1*z3-
q1*z4)*QDSx([3, 5, 2, 1, 4], 'y') + q1*QDSx([3, 5, 4, 1, 2], 'y') + (q1*y1+q1*y3-q1*z1-q1*z2)*QDSx([3, 
6, 1, 2, 4, 5], 'y') + q1*QDSx([3, 6, 2, 1, 4, 5], 'y') + (q3*y4^4-q3*y4^3*z1-q3*y4^3*z2+q3*y4^2*z1*z2-
q3*y4^3*z3+q3*y4^2*z1*z3+q3*y4^2*z2*z3-q3*y4*z1*z2*z3-q3*y4^3*z4+q3*y4^2*z1*z4+q3*y4^2*z2*z4-
q3*y4*z1*z2*z4+q3*y4^2*z3*z4-q3*y4*z1*z3*z4-q3*y4*z2*z3*z4+q3*z1*z2*z3*z4)*QDSx([4, 1, 2, 3], 'y') + 
(y1*y4^4+y3*y4^4-y1*y4^3*z1-y3*y4^3*z1-y4^4*z1+y4^3*z1^2-y1*y4^3*z2-y3*y4^3*z2-y4^4*z2+y1*y4^2*z1*z2+
y3*y4^2*z1*z2+2*y4^3*z1*z2-y4^2*z1^2*z2+y4^3*z2^2-y4^2*z1*z2^2-y1*y4^3*z3-y3*y4^3*z3+y1*y4^2*z1*z3+
y3*y4^2*z1*z3+y4^3*z1*z3-y4^2*z1^2*z3+y1*y4^2*z2*z3+y3*y4^2*z2*z3+y4^3*z2*z3-y1*y4*z1*z2*z3-
y3*y4*z1*z2*z3-2*y4^2*z1*z2*z3+y4*z1^2*z2*z3-y4^2*z2^2*z3+y4*z1*z2^2*z3-y1*y4^3*z4-y3*y4^3*z4+
y1*y4^2*z1*z4+y3*y4^2*z1*z4+y4^3*z1*z4-y4^2*z1^2*z4+y1*y4^2*z2*z4+y3*y4^2*z2*z4+y4^3*z2*z4-
y1*y4*z1*z2*z4-y3*y4*z1*z2*z4-2*y4^2*z1*z2*z4+y4*z1^2*z2*z4-y4^2*z2^2*z4+y4*z1*z2^2*z4+y1*y4^2*z3*z4+
y3*y4^2*z3*z4-y1*y4*z1*z3*z4-y3*y4*z1*z3*z4-y4^2*z1*z3*z4+y4*z1^2*z3*z4-y1*y4*z2*z3*z4-y3*y4*z2*z3*z4-
y4^2*z2*z3*z4+y1*z1*z2*z3*z4+y3*z1*z2*z3*z4+2*y4*z1*z2*z3*z4-z1^2*z2*z3*z4+y4*z2^2*z3*z4-
z1*z2^2*z3*z4)*QDSx([4, 1, 3, 2], 'y') + (y4^4-y4^3*z1-y4^3*z2+y4^2*z1*z2-y4^3*z3+y4^2*z1*z3+y4^2*z2*z3-
y4*z1*z2*z3-y4^3*z4+y4^2*z1*z4+y4^2*z2*z4-y4*z1*z2*z4+y4^2*z3*z4-y4*z1*z3*z4-y4*z2*z3*z4+
z1*z2*z3*z4)*QDSx([4, 1, 5, 2, 3], 'y') + (y4^4-y4^3*z1-y4^3*z2+y4^2*z1*z2-y4^3*z3+y4^2*z1*z3+y4^2*z2*z3-
y4*z1*z2*z3-y4^3*z4+y4^2*z1*z4+y4^2*z2*z4-y4*z1*z2*z4+y4^2*z3*z4-y4*z1*z3*z4-y4*z2*z3*z4+
z1*z2*z3*z4)*QDSx([4, 2, 3, 1], 'y') + (q1*y1+q1*y3+q1*y4+q1*y5-q1*z1-q1*z2-q1*z3-q1*z4)*QDSx([4, 5, 1, 
2, 3], 'y') + q1*QDSx([4, 5, 2, 1, 3], 'y') + q1*QDSx([4, 6, 1, 2, 3, 5], 'y') + (q3*y4^3+q3*y4^2*y5+
q3*y4*y5^2+q3*y5^3-q3*y4^2*z1-q3*y4*y5*z1-q3*y5^2*z1-q3*y4^2*z2-q3*y4*y5*z2-q3*y5^2*z2+q3*y4*z1*z2+
q3*y5*z1*z2-q3*y4^2*z3-q3*y4*y5*z3-q3*y5^2*z3+q3*y4*z1*z3+q3*y5*z1*z3+q3*y4*z2*z3+q3*y5*z2*z3-
q3*z1*z2*z3-q3*y4^2*z4-q3*y4*y5*z4-q3*y5^2*z4+q3*y4*z1*z4+q3*y5*z1*z4+q3*y4*z2*z4+q3*y5*z2*z4-
q3*z1*z2*z4+q3*y4*z3*z4+q3*y5*z3*z4-q3*z1*z3*z4-q3*z2*z3*z4)*QDSx([5, 1, 2, 3, 4], 'y') + (y1*y4^3+
y3*y4^3+y1*y4^2*y5+y3*y4^2*y5+y1*y4*y5^2+y3*y4*y5^2+y1*y5^3+y3*y5^3-y1*y4^2*z1-y3*y4^2*z1-y4^3*z1-
y1*y4*y5*z1-y3*y4*y5*z1-y4^2*y5*z1-y1*y5^2*z1-y3*y5^2*z1-y4*y5^2*z1-y5^3*z1+y4^2*z1^2+y4*y5*z1^2+
y5^2*z1^2-y1*y4^2*z2-y3*y4^2*z2-y4^3*z2-y1*y4*y5*z2-y3*y4*y5*z2-y4^2*y5*z2-y1*y5^2*z2-y3*y5^2*z2-
y4*y5^2*z2-y5^3*z2+y1*y4*z1*z2+y3*y4*z1*z2+2*y4^2*z1*z2+y1*y5*z1*z2+y3*y5*z1*z2+2*y4*y5*z1*z2+
2*y5^2*z1*z2-y4*z1^2*z2-y5*z1^2*z2+y4^2*z2^2+y4*y5*z2^2+y5^2*z2^2-y4*z1*z2^2-y5*z1*z2^2-y1*y4^2*z3-
y3*y4^2*z3-y1*y4*y5*z3-y3*y4*y5*z3-y1*y5^2*z3-y3*y5^2*z3+y1*y4*z1*z3+y3*y4*z1*z3+y4^2*z1*z3+y1*y5*z1*z3+
y3*y5*z1*z3+y4*y5*z1*z3+y5^2*z1*z3-y4*z1^2*z3-y5*z1^2*z3+y1*y4*z2*z3+y3*y4*z2*z3+y4^2*z2*z3+y1*y5*z2*z3+
y3*y5*z2*z3+y4*y5*z2*z3+y5^2*z2*z3-y1*z1*z2*z3-y3*z1*z2*z3-2*y4*z1*z2*z3-2*y5*z1*z2*z3+z1^2*z2*z3-
y4*z2^2*z3-y5*z2^2*z3+z1*z2^2*z3-y1*y4^2*z4-y3*y4^2*z4-y1*y4*y5*z4-y3*y4*y5*z4-y1*y5^2*z4-y3*y5^2*z4+
y1*y4*z1*z4+y3*y4*z1*z4+y4^2*z1*z4+y1*y5*z1*z4+y3*y5*z1*z4+y4*y5*z1*z4+y5^2*z1*z4-y4*z1^2*z4-y5*z1^2*z4+
y1*y4*z2*z4+y3*y4*z2*z4+y4^2*z2*z4+y1*y5*z2*z4+y3*y5*z2*z4+y4*y5*z2*z4+y5^2*z2*z4-y1*z1*z2*z4-
y3*z1*z2*z4-2*y4*z1*z2*z4-2*y5*z1*z2*z4+z1^2*z2*z4-y4*z2^2*z4-y5*z2^2*z4+z1*z2^2*z4+y1*y4*z3*z4+
y3*y4*z3*z4+y1*y5*z3*z4+y3*y5*z3*z4-y1*z1*z3*z4-y3*z1*z3*z4-y4*z1*z3*z4-y5*z1*z3*z4+z1^2*z3*z4-
y1*z2*z3*z4-y3*z2*z3*z4-y4*z2*z3*z4-y5*z2*z3*z4+2*z1*z2*z3*z4+z2^2*z3*z4)*QDSx([5, 1, 3, 2, 4], 'y') + 
(y4^3+y4^2*y5+y4*y5^2+y5^3-y4^2*z1-y4*y5*z1-y5^2*z1-y4^2*z2-y4*y5*z2-y5^2*z2+y4*z1*z2+y5*z1*z2-y4^2*z3-
y4*y5*z3-y5^2*z3+y4*z1*z3+y5*z1*z3+y4*z2*z3+y5*z2*z3-z1*z2*z3-y4^2*z4-y4*y5*z4-y5^2*z4+y4*z1*z4+y5*z1*z4+
y4*z2*z4+y5*z2*z4-z1*z2*z4+y4*z3*z4+y5*z3*z4-z1*z3*z4-z2*z3*z4)*QDSx([5, 1, 4, 2, 3], 'y') + (y4^3+
y4^2*y5+y4*y5^2+y5^3-y4^2*z1-y4*y5*z1-y5^2*z1-y4^2*z2-y4*y5*z2-y5^2*z2+y4*z1*z2+y5*z1*z2-y4^2*z3-
y4*y5*z3-y5^2*z3+y4*z1*z3+y5*z1*z3+y4*z2*z3+y5*z2*z3-z1*z2*z3-y4^2*z4-y4*y5*z4-y5^2*z4+y4*z1*z4+y5*z1*z4+
y4*z2*z4+y5*z2*z4-z1*z2*z4+y4*z3*z4+y5*z3*z4-z1*z3*z4-z2*z3*z4)*QDSx([5, 2, 3, 1, 4], 'y') + (q3*y4^2+
q3*y4*y5+q3*y5^2+q3*y4*y6+q3*y5*y6+q3*y6^2-q3*y4*z1-q3*y5*z1-q3*y6*z1-q3*y4*z2-q3*y5*z2-q3*y6*z2+
q3*z1*z2-q3*y4*z3-q3*y5*z3-q3*y6*z3+q3*z1*z3+q3*z2*z3-q3*y4*z4-q3*y5*z4-q3*y6*z4+q3*z1*z4+q3*z2*z4+
q3*z3*z4)*QDSx([6, 1, 2, 3, 4, 5], 'y') + (y1*y4^2+y3*y4^2+y1*y4*y5+y3*y4*y5+y1*y5^2+y3*y5^2+y1*y4*y6+
y3*y4*y6+y1*y5*y6+y3*y5*y6+y1*y6^2+y3*y6^2-y1*y4*z1-y3*y4*z1-y4^2*z1-y1*y5*z1-y3*y5*z1-y4*y5*z1-y5^2*z1-
y1*y6*z1-y3*y6*z1-y4*y6*z1-y5*y6*z1-y6^2*z1+y4*z1^2+y5*z1^2+y6*z1^2-y1*y4*z2-y3*y4*z2-y4^2*z2-y1*y5*z2-
y3*y5*z2-y4*y5*z2-y5^2*z2-y1*y6*z2-y3*y6*z2-y4*y6*z2-y5*y6*z2-y6^2*z2+y1*z1*z2+y3*z1*z2+2*y4*z1*z2+
2*y5*z1*z2+2*y6*z1*z2-z1^2*z2+y4*z2^2+y5*z2^2+y6*z2^2-z1*z2^2-y1*y4*z3-y3*y4*z3-y1*y5*z3-y3*y5*z3-
y1*y6*z3-y3*y6*z3+y1*z1*z3+y3*z1*z3+y4*z1*z3+y5*z1*z3+y6*z1*z3-z1^2*z3+y1*z2*z3+y3*z2*z3+y4*z2*z3+
y5*z2*z3+y6*z2*z3-2*z1*z2*z3-z2^2*z3-y1*y4*z4-y3*y4*z4-y1*y5*z4-y3*y5*z4-y1*y6*z4-y3*y6*z4+y1*z1*z4+
y3*z1*z4+y4*z1*z4+y5*z1*z4+y6*z1*z4-z1^2*z4+y1*z2*z4+y3*z2*z4+y4*z2*z4+y5*z2*z4+y6*z2*z4-2*z1*z2*z4-
z2^2*z4+y1*z3*z4+y3*z3*z4-z1*z3*z4-z2*z3*z4)*QDSx([6, 1, 3, 2, 4, 5], 'y') + (y4^2+y4*y5+y5^2+y4*y6+
y5*y6+y6^2-y4*z1-y5*z1-y6*z1-y4*z2-y5*z2-y6*z2+z1*z2-y4*z3-y5*z3-y6*z3+z1*z3+z2*z3-y4*z4-y5*z4-y6*z4+
z1*z4+z2*z4+z3*z4)*QDSx([6, 1, 4, 2, 3, 5], 'y') + (y4^2+y4*y5+y5^2+y4*y6+y5*y6+y6^2-y4*z1-y5*z1-y6*z1-
y4*z2-y5*z2-y6*z2+z1*z2-y4*z3-y5*z3-y6*z3+z1*z3+z2*z3-y4*z4-y5*z4-y6*z4+z1*z4+z2*z4+z3*z4)*QDSx([6, 2, 
3, 1, 4, 5], 'y') + (q3*y4+q3*y5+q3*y6+q3*y7-q3*z1-q3*z2-q3*z3-q3*z4)*QDSx([7, 1, 2, 3, 4, 5, 6], 'y') + 
(y1*y4+y3*y4+y1*y5+y3*y5+y1*y6+y3*y6+y1*y7+y3*y7-y1*z1-y3*z1-y4*z1-y5*z1-y6*z1-y7*z1+z1^2-y1*z2-y3*z2-
y4*z2-y5*z2-y6*z2-y7*z2+2*z1*z2+z2^2-y1*z3-y3*z3+z1*z3+z2*z3-y1*z4-y3*z4+z1*z4+z2*z4)*QDSx([7, 1, 3, 2, 
4, 5, 6], 'y') + (y4+y5+y6+y7-z1-z2-z3-z4)*QDSx([7, 1, 4, 2, 3, 5, 6], 'y') + (y4+y5+y6+y7-z1-z2-z3-
z4)*QDSx([7, 2, 3, 1, 4, 5, 6], 'y') + q3*QDSx([8, 1, 2, 3, 4, 5, 6, 7], 'y') + (y1+y3-z1-z2)*QDSx([8, 
1, 3, 2, 4, 5, 6, 7], 'y') + QDSx([8, 1, 4, 2, 3, 5, 6, 7], 'y') + QDSx([8, 2, 3, 1, 4, 5, 6, 7], 
'y')
```
This output of roughly 17,000 characters took about 2 seconds to compute on my laptop. **Again note that 
quantum double computations are technically conjectural, but a proof is likely forthcoming.**


[Homepage of schubmult](http://schubmult.org/)