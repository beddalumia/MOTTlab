## MIT tutorial

GOAL: to run the IPT codes for the (approximate) solution of 
      the Hubbard model in DMFT


* code:  pertre for real frequency IPT solution


* `input' file: fort.22 for pertre.f

download the 2 files in your account.


The files in detail:

***********************************************************

Input file for real axis IPT:

fort.22

0.0005,1.,2.3

dE, D, U

dE energy mesh
D half-bandwidth
U Hubbard interaction
* * * * * * * * * * 


************************************************************

output files:

fort.25 is Im[G(w)]   (the DOS up to a factor of pi)
fort.23 is Re[Sigma(w)]
fort.24 is Im[Sigma(w)]
fort.26 is Im[G(w)]  (impurity problem, no self-consistency, imet=1)

************************************************************


************************************************************


Lets run the codes now

Excercises

1. Run the real-frequency code pertre.f for several values of U (set the 
half-bandwidth D=1) to check out the metal-insulator transition at T=0.
The DOS(w) is written in "fort.25", the Selfenergy in "fort.23" and "fort.24"
(real and imaginary parts). Note that Re[Sigma] is linear for a metal (Fermi liquid)
and that it diverges as 1/w in the insulator.

2. Approach the MIT at Uc_2 \approx 3.D (metal to insulator) from below and observe 
how the different contributions to the DOS (quasiparticle peak and Hubbard
bands) evolve as you get close to the critical value.

3. Check for the existence of 2 different solutions
When you are within the coexistence region, at Uc_1 < U <Uc_2 
(Uc_1 \approx 2.6, Uc_2 \approx 3.) the code will ask you to choose either the metallic 
or the insulating solutions. You can solve for one and save the output and the solve for 
the other solution and compare the solutions.

3. The effect of the self-consistency condition: For several values of
U (U < Uc_2), run the real-frequency code and compare the output files
"fort.25" and "fort.26" that contain the converged DOS(w) and the DOS(w)
after only the first iteration. 
Thus, "fort.26" corresponds to the solution of the Single Impurity Anderson Model 
for a semi-circular conduction band of total bandwidth W=2D=2.

Comparing the outputs fort.25 and fort.26, you can thus observe the effect of the 
self-consistency. 
What is the main effect?
Thinking in terms of the SIAM (or Kondo model) can you explain qualitatively the 
behavior of the quasiparticle peak of the Hubbard model?


