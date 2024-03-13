# Overview

This is the implementation of a deterministic finite automata that generates a password of length at least *k*, accepted by the automata.

The automata has *n* states, numbered from 1 to *n* (*q<sub>1</sub>, q<sub>2</sub>, ..., q<sub>n</sub>*). Also, the entry alphabet is $\Sigma$, the initial state is *q<sub>x</susb>*, and there are *m* final states *{q<sub>s<sub>1</sub></sub>, q<sub>s<sub>2</sub></sub>, ..., q<sub>s<sub>m</sub></sub>}* with the transition state $\delta$.

The scope of this automata is to find a password *w* that has the least length possible over *k*.

It is implied that $\Sigma$ is formed out of the first |$\Sigma$| lowercased letter of the latin alphabet.

## Input

The input data is written in a file named **input.txt** having the following format.

The first line contains a number *t* $\epsilon$ {1, 2} which specifies the subtask run by the program.
The difference between the two is based on the size of the problem to be solved. 
The first subtask uses a simpler algorithm, while the second uses a more complex one in order to solve more complex problems. 
For more details look at the constraints.

The second line contains 4 numbers: *n*, *|$\Sigma$|*, *m*, and *k*.
Their significance is described above.

The third line contains a single number, *q<sub>x</sub>*, the initial state.

The forth line contains *m* numbers, *{q<sub>s<sub>1</sub></sub>, q<sub>s<sub>2</sub></sub>, ..., q<sub>s<sub>m</sub></sub>}*, the final states.

The last *n* lines contain the rules for the $\delta$ function.

## Output

The output is viewed in **output.txt**.
If there is no solution found, it will be displayed -1.

Otherwise, on the first line, it should be displayed *|w|*, the lenght of the password, followed by the password on the second line.

## Constraints

* $2 \leq |\Sigma| \leq 26$, $2 \leq n$, $1 \leq k$
* $1 \leq$ *q<sub>x</sub>* $\leq n.$ For any *i, j*, $1 \leq i < j \leq m,$ *q<sub>s<sub>i</sub></sub>* $\neq$ *q<sub>s<sub>j</sub></sub>*  
* Solutions are accepted only for *|w|* $\leq 2 \cdot 10^8$

Subproblem 1:

* $n \leq 5000, k \leq 2 \cdot 10^4$

Subproblem 2:

* $n \leq 200, k \leq 10^8$
