### Introduction

This is the 4th task of CMM.

In this task, I implement a 2D FEM simulator based on templates from TAs.

More details, check "tutorials-a4" and "FESA".

### Demo

check them in "/media"

![fem](https://user-images.githubusercontent.com/39910677/116280919-76c65680-a789-11eb-86a1-dc9a769d5ecc.gif) ![manip](https://user-images.githubusercontent.com/39910677/116280986-880f6300-a789-11eb-8c3c-a229c56db4c5.gif)






### How to start

```git clone```

```cd 2D_FEM_Simulation/build/src/app```

```./fem-app```

```./manip-app```

### Solution to Question 10:
Consider the loss function: $O(u) = 0.5*||x(u)-x'||^2 : R^6 \rightarrow R $ For a given fixed pin position pair x', there exists at least two different u vectors (shown in the figure) s.t. O(u) = 0. And because O(u) is non-negative, the two u are global minima. So the function is non-convex.

The algorithm to find a solution is gradient descent, which is sensitive to the initial value. So the optimizer tends to find the first solution because of the initial u configuration.

### Solution to Question 11:
To solve this problem, we should apply a regularizer to the optimized loss function. The intuitive idea is to add a norm punishment to vector u, but it performs poorly. Then I try to add the total bar energy as the regularizer.

After several trials, I find that the stress in the first solution is not uniform. So I turn to use the difference of maximum stress and minimum stress as the regularizer. It works better than total energy.

Please check the code for implementation details.
