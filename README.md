# Vicsek model

C++ code on numerical simulation of the Vicsek model. The algorithm used here is explained in the article: F. Ginelli, <a href='https://link.springer.com/article/10.1140/epjst/e2016-60066-8'><i>The Physics of the Vicsek model</i></a>, Eur. Phys. J. Spec. Top. <b>225</b>, 2099–2117 (2016). The model was studied in refs:<br>
[1] T. Vicsek, A. Czirók, E. Ben-Jacob, I. Cohen, and Ofer Shochet, <a href='https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.75.1226'><i>Novel Type of Phase Transition in a System of Self-Driven Particles</i></a>, Phys. Rev. Lett. <b>75</b>, 1226 (1995);<br>
[2] A. P. Solon, H. Chaté, and J. Tailleur, <a href='https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.068101'><i>From Phase to Microphase Separation in Flocking Models: The Essential Role of Nonequilibrium Fluctuations</i></a>, Phys. Rev. Lett. <b>114</b>, 068101 (2015);<br>
[3] A. P. Solon, J.-B. Caussin, D. Bartolo, Hugues Chaté, and J. Tailleur, <a href='https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.062111'><i>Pattern formation in flocking models: A hydrodynamic description</i></a>, Phys. Rev. E 92, 062111 (2015);<br>
[4] R. Kürsten and T. Ihle, <a href='https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.188003'><i>Dry Active Matter Exhibits a Self-Organized Cross Sea Phase</i></a>, Phys. Rev. Lett. <b>125</b>, 188003 (2020).<br>
<br>
Exportations: dynamics of flocking; time evolution of the order parameters and the mean-square displacement.<br>
Compile: g++ VM.cpp -lgsl -lgslcblas -lm -O3 -s -o VM.out.<br>
Run: ./VM.out -parameter=value.<br>
List of parameters: v0, eta, rho0, LX, LY, tmax, init, RAN (details as comments in the code).<br>
Generate the movie (rectangle geometry 800x100): python figure_VM_dynamics.py -parameter=value (parameters: v0, eta, rho0, LX, LY, DT, tmax, init, ran).<br>
Generate the movie (square geometry 512x512): python figure_VM_square_dynamics.py -parameter=value (parameters: v0, eta, rho0, LX, LY, DT, tmax, init, ran).<br>
