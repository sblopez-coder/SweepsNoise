# Created by Ben Blowers as a part of the Fermilab SQMS summer internship, 2021

from qutip import *
import numpy as np
import scqubits as scq
from matplotlib import pyplot, animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import parameters

scq.settings.T1_DEFAULT_WARNING=False # suppresses some warning

# --Set Parameters--
Ej = parameters.Ej
Ec = parameters.Ec
ng = parameters.ng
cavity = parameters.cavity
atom = parameters.atom
g = parameters.g
# --/Set Parameters--

sample = 400

# --Take User Inputs
print("What sweep to do?")
print("{1} Energy States")
print("{2} Transition Energies")
print("{3} Effective T1")
print("{4} Rabi Oscillations")
type = int(input())
if (type != 1 and type != 2 and type != 3 and type != 4):
    print("Yo, get your act together. " + str(type) + " wasn't an option.")
    exit()

if (type == 1 or type == 2 or type == 3):
    print()
    print("Plot against which parameter:")
    print("1: Offset Charge (ng)")
    print("2: Josephson Energy (Ej)")
    print("3: Charging Energy (Ec)")
    parameter = int(input())
    if (parameter != 1 and parameter != 2 and parameter != 3):
        print("Yo, get your act together. " + str(parameter) + " wasn't an option.")
        exit()
    elif (parameter == 1):
        parameter = "ng"
    elif (parameter == 2):
        parameter = "EJ"
    elif (parameter == 3):
        parameter = "EC"
    print()
    range = int(input("Which range of parameter values to sweep over?"))
    if (type == 1 or type == 3):
        params = np.linspace(-range,range,sample)
    if (type == 2):
        params = np.linspace(0,range,sample)

if (type == 4):
    mode_p = int(input("Vary the {1} amplitude, {2} length, or {3} both?"))
    if (mode_p == 1 or mode_p == 2):
        range = int(input("Enter the range of values to sweep over: "))
        if (mode_p == 1):
            params = np.linspace(-range, range, sample)
            width = float(input("Enter the pulse width (this will be constant): "))
        if (mode_p == 2):
            params = np.linspace(0, range, sample)
            A = float(input("Enter the pulse amplitude (this will be constant): "))
        print()
        makeSphere = int(input("Make a Bloch sphere animation? {0} for no, {1} for yes"))
    if (mode_p == 3):
        rangeA = int(input("Enter the range of amplitudes to sweep over: "))
        amplitudes = np.linspace(-rangeA,rangeA, sample)
        rangeW = int(input("Enter the range of lengths to sweep over: "))
        widths = np.linspace(0,rangeW, sample)
        makeSphere = 0 # Making sure that this is declared

# Create the qubit
qubit = scq.Transmon(EJ=Ej, EC=Ec, ng=ng, ncut=150) # ncut seems to be just a number that needs to be big
N = 2
wc = cavity*2*np.pi
wq = atom*2*np.pi
a = tensor(qeye(2),destroy(N))
sm = tensor(destroy(2),qeye(N))
sx = tensor(sigmax(),qeye(N))
sy = tensor(sigmay(),qeye(N))
sz = tensor(sigmaz(),qeye(N))
ground = tensor(basis(2,0), basis(N,0))
excite = tensor(basis(2,1), basis(N,0))

# Defining the drive pulses
def square(t, args):
    return args['A']*(t > args['delay']) - args['A']*(t>(args['width']+args['delay']))

def gauss(t, args):
    return args['A'] * np.exp(-((t-.5) / args['width']) ** 2)

# Hamiltonians
H0 = -0.5*wq*sz # qubit Hamiltonian
Hc = wc*(a.dag()*a) # cavity Hamiltonian
Hi = g*(a.dag()*sm + a*sm.dag()) # interaction Hamiltonian
Hdrive = a+a.dag() # drive Hamiltonian
H = [H0,[Hdrive,square]] # replace 'square' with 'gauss' to use a gaussian pulse

c_ops = []
e_ops = [excite*excite.dag()]

# this Rabi sweep works by running mesolve over and over again, using a different amplitude and/or length each time
# the resulting states are stored into a 'results' array
if (type == 4):
    if (mode_p == 1 or mode_p == 2):
        results = np.zeros(len(params))
        for i in np.arange(0,len(params)):
            if (mode_p == 1):
                tlist = np.linspace(0, width + 2, 10)
                result = mesolve(H, excite, tlist, c_ops, e_ops, args={'A': params[i], 'width': width, 'delay': 0}, options=Options(nsteps=5000))
            if (mode_p == 2):
                tlist = np.linspace(0, params[i] + 2, 10)
                result = mesolve(H, excite, tlist, c_ops, e_ops, args={'A': A, 'width': params[i], 'delay': 0}, options=Options(nsteps=5000))
            results[i] = 1 - result.expect[0][len(tlist)-1]
            print(results[i])

        fig, ax = plt.subplots(figsize=(12, 6))
        plt.rcParams.update({'font.size': 22})
        ax.plot(params, results, 'b')
        if (mode_p == 1):
            ax.set_xlabel('Amplitude', fontsize=20)
        if (mode_p == 2):
            ax.set_xlabel('Length', fontsize=20)
        ax.set_ylabel('Probability', fontsize=20)
        plt.ylim([-.1, 1.1])
        plt.show()

    if (mode_p == 3):
        tlist = np.linspace(0, widths[sample - 1] + 5, 50)
        results = [[0 for x in np.arange(0,len(amplitudes))] for y in np.arange(0,len(widths))]
        for i in np.arange(0,len(amplitudes)):
            for j in np.arange(0,len(widths)):
                result = mesolve(H, excite, tlist, c_ops, e_ops, args={'A': amplitudes[i], 'width': widths[j], 'delay': 0},options=Options(nsteps=5000))
                results[j][i] = result.expect[0][len(tlist) - 1]
                print(str(i) + ", " + str(j) + ": " + str(results[i][j]))

        # rotating the axes for the 2D plot
        Results = [[0 for x in np.arange(0,len(amplitudes))] for y in np.arange(0,len(widths))]
        for i in np.arange(0,len(results)):
            for j in np.arange(0,len(results[0])):
                Results[i][j] = 1-results[sample-1-i][j]

        fig, ax = plt.subplots(figsize=(12,6))
        img = ax.imshow(Results,cmap = 'plasma', interpolation='none')
        plt.title('Rabi Sweeps')
        ax.set_xlabel('Amplitude')
        ax.set_ylabel('Length')
        plt.xticks([sample/4, sample/2, sample*3/4],[-rangeA/2, 0, rangeA/2])
        plt.yticks([sample/4, sample/2, sample*3/4],[rangeW*3/4, rangeW/2, rangeW/4])
        fig.colorbar(img)
        plt.savefig("rabi sweeps.png")
        plt.show()

    if (makeSphere == 1):
        fig = pyplot.figure()
        ax = Axes3D(fig, azim=-40, elev=30)
        sphere = qutip.Bloch(axes=ax)

        # Convert expectation values to spherical coordinates
        phase = np.diff(results)
        phase = np.append(phase, phase[len(tlist)-1])
        theta = results * np.pi
        phi = np.zeros(len(phase))
        for i in np.arange(len(phase)):
            if phase[i] >= 0:
                phi[i] = 0
            else:
                phi[i] = np.pi
        print(phi)

        def animate(i):
            sphere.clear()
            sphere.add_vectors([np.sin(theta[i]) * np.cos(phi[i]), np.sin(theta[i]) * np.sin(phi[i]), np.cos(theta[i])])
            sphere.make_sphere()
            return ax

        def init():
            sphere.vector_color = ['r']
            return ax

        ani = animation.FuncAnimation(fig, animate, np.arange(len(results)), init_func=init, blit=False, repeat=False)
        ani.save("rabi.gif", fps=15)

if (type == 1):
    qubit.plot_evals_vs_paramvals(parameter, params, evals_count=5, subtract_ground=False)
    plt.show()

if (type == 2):
    qubit.plot_dispersion_vs_paramvals('ng', parameter, params, transitions=(((0, 1), (0, 2))))
    plt.show()

if (type == 3):
    qubit.plot_t1_effective_vs_paramvals(parameter, params)
    plt.show()