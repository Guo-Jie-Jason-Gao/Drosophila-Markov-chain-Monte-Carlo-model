# Drosophila-Markov-chain-Monte-Carlo-model
Code for our manuscript "A Markov chain Monte Carlo model of
mechanical-feedback-driven progressive
apical constrictions captures the fluctuating
collective cell dynamics in the Drosophila embryo" 
by Guo–Jie J. Gao, Michael C. Holcomb, Jeffrey H. Thomas, and Jerzy Blawzdziewicz.

To run the code:
Step 1) Compile main.f by typing "ifort main.f" on a unix / linux / Mac computer.
An Intel Fortran compiler is needed.

Step 2) Run the executable by typing 
"./a.out < input_file". The executable will read an initial configuration "IC" automatically.
