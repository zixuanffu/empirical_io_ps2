# Empirical Industrial Organization Problem Set 2
## Dynamic Discrete Choice 
In this assignment you have to simulate a simple version of an inventory control problem. The objective
is to familiarize with the computation of the dynamic programming solution and compare it with a CCP
method.

Deadline: Please send your problem set before January 10th, 2025.
## Requirements
see [here](/Reports/PS2_corrected.pdf)  

## Remarks
There are two types of penalty term in the utility function, each with their economic interpretation.
### Case 1: $-\lambda 1\{c>0\} 1\{i=0\}$
1. In this period t, I somehow want to eat chips $c>0$, but I have run out of chips $i=0$.
2. **I incur some penalty because I can not eat my fav snack right away.** I need to go out and buy some.
3. Then I go to the supermarket, see the price $(p_s)$. I know the state $(i=0,c>0,p, \epsilon)$ when making the decision $x$. 
4. If I decide to buy $x=1$, then my next period inventory would be $i'=i+x-c=0.75$
5. If I decide not to buy $x=0$, then my next period $i'=\max\{0,i+x-c\}=0$

### Case 2: $-\lambda 1\{c>0\} 1\{i+x=0\}$
1. In this period t, I somehow want to eat chips $c>0$, but I have run out of chips $i=0$.
2. Then I go to the supermarket, see the price $(p_s)$. I know the state $(i=0,c>0,p, \epsilon)$ when making the decision $x$. 
3. If I decide to buy $x=1$, then I can happily eat the chips and there's no penalty. My next period $i'=i+x-c=0.75$
4. If I decide not to buy $x=0$, I incur the loss of "starving". Then my next period $i'=\max\{0,i+x-c\}=0$
## Code
see [here](/Code/)
## Report
see [here](/Reports/report.pdf)
