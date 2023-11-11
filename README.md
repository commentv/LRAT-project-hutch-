# Presentation
This project aims on comparing two trace estimators of a large symmetric matrix A which we can only access via matrix-vector multiplications.
The first method to estimate $tr(A)$ is using a Monte-Carlo estimator, where, given $m$ i.i.d. gaussian vectors $w_1, ... , w_m \sim N(0,I_n)$,
then :<br />
<p align="center">
$$tr_m(A) = \sum_{i=1}^m w_i^tAw_i $$
</p>
The second method is inspired from [1, Algorithm 1] but with random Gaussian vectors instead of Rademacher vectors. The codes and the results follow the questions
given in the ```Hutch++_questions.pdf``` file.

# Explanation
To generate the plots used to answer questions 1, 3 and 4 of the ```Hutch++_questions.pdf``` file and stored in ```/plots```, run their respective matlab code: 
```/src/Question1.m```, ```/src/Question3.m``` and ```/src/Question4.m```.<br /><br />
Note that the plots were generated using Matlab R2022b and the library "Statistics and Machine Learning Toolbox" is needed.
Moreover, the plot from ```/src/Question4.m``` answers question 4 by giving an example on how the decay might be slower than $O(\frac{1}{m})$ while respecting the inequalities
of Theorem 1 and Theorem 10 of [1].

# References
[1] Meyer, R., Musco, C., Musco, C., Woodruff, D.P. Hutch++: Optimal
Stochastic Trace Estimation arXiv preprint arXiv:2010.09649 (2021). https:
//arxiv.org/pdf/2010.09649
