import numpy as np
import matplotlib.pyplot as plt 

initialVelocity = 10 # m/s
initialAngle = 30 # degrees
time = 40 # s
tstep = 0.01 # time step

# M80 FMJ Ball
r = 0.00381 # m Bullet radius (7.62x51 mm NATO) 0.3085in
ogiveLength = 0.0194818 # m 0.767in
ogiveRadius = 0.0762 # m 3.00in
meplatDiameter = 0.001524 # m 0.060in
boatTailLength = 0.004445 # m .175in

barrelTwistRate = 12 # Inches/rotation

boatTailAngle = np.radians(9.0) # deg
oal = 0.028956 # m 1.140in
bulletMass = 9.52544e-3 # kg 9.52544g = 147gr

# Assume homogenous density and approximate mass distribution with cones and cylinders
areaSection = np.pi*r**2
BaseRadius = r - boatTailLength/np.tan(boatTailAngle)
volumeBaseConeSmall = np.pi*BaseRadius**3*np.tan(boatTailAngle)/3
volumeBaseConeBig = np.pi*r**3*np.tan(boatTailAngle)/3
volumeOgive = np.pi*r**2*ogiveLength/3
volumeBoatTail = (np.pi*np.tan(boatTailAngle)/3)*(r**3 - BaseRadius**3)
volumeBearing = np.pi*r**2*(oal - boatTailLength - ogiveLength)
volumeTotal = volumeOgive + volumeBoatTail + volumeBearing
density = bulletMass/volumeTotal
massOgive = volumeOgive * density
massBoatTail = volumeBoatTail * density
massBearing = volumeBearing * density
massBaseConeSmall = volumeBaseConeSmall*density
massBaseConeBig = volumeBaseConeBig*density
IogiveInline = 3*massOgive*r**2/10
IbtConeSmallInline = 3*massBaseConeSmall*BaseRadius**2/10
IbtConeBigInline = 3*massBaseConeBig*r**2/10
IboatTailInline = IbtConeBigInline - IbtConeSmallInline
IbearingInline = massBearing*r**2/2
Inertia_Inline = IogiveInline + IboatTailInline + IbearingInline
IbearingOffAxis = massBearing*(r**2/4 + ((oal - ogiveLength - boatTailLength)**2)/12)
IbtConeSmallOffAxis = massBaseConeSmall*((BaseRadius*np.tan(boatTailAngle))**2/10 + 3*BaseRadius**2/20)
IbtConeBigOffAxis = massBaseConeBig*((r*np.tan(boatTailAngle))**2/10 + 3*r**2/20)
InertiaBoatTailOffAxis = IbtConeBigOffAxis - IbtConeSmallOffAxis
InertiaOgiveOffAxis = massOgive*(ogiveLength**2/10 + 3*r**2/20)
Inertia_OffAxis = IbearingOffAxis + InertiaBoatTailOffAxis + InertiaOgiveOffAxis

initialSpin = (12/barrelTwistRate)*(initialVelocity*(2*np.pi*3.28084))


def ode45py(func, x, y, st_sz=1.0e-4, tol=1.0e-6, iter_lim=50000):
    '''
    Numerical Methods: Differential Equations, Initial Value Problems

    4th-order / 5th-order Runge-Kutta Method
    Includes adaptive step size adjustment
    Imitates MATLAB ode45 functionality and output
    '''
    # Dormand-Prince coefficients for RK algorithm -
    a1 = 0.2; a2 = 0.3; a3 = 0.8; a4 = 8/9; a5 = 1.0; a6 = 1.0
    c0 = 35/384; c2 = 500/1113; c3 = 125/192; c4 = -2187/6784; c5=11/84
    d0 = 5179/57600; d2 = 7571/16695; d3 = 393/640; d4 = -92097/339200; d5 = 187/2100; d6 = 1/40
    b10 = 0.2
    b20 = 0.075; b21 = 0.225
    b30 = 44/45; b31 = -56/15; b32 = 32/9
    b40 = 19372/6561; b41 = -25360/2187; b42 = 64448/6561; b43 = -212/729
    b50 = 9017/3168; b51 = -355/33; b52 = 46732/5247; b53 = 49/176; b54 = -5103/18656
    b60 = 35/384; b62 = 500/1113; b63 = 125/192; b64 = -2187/6784; b65 = 11/84
    # Store initial values
    x_f = x[-1]
    x_n = x[0]
    # y_n = y
    # Initialize variables
    X = []
    Y = []
    # Add the first set of known conditions
    X.append(x_n)
    Y.append(y)
    # Set up to break the for loop at the end
    stopper = 0 # Integration stopper, 0 = off, 1 = on
    # Initialize a k0 to start with the step size
    k0 = st_sz * func(x_n, y)
    # Generate the RK coefficients
    for i in range(iter_lim):
        k1 = st_sz * func(x_n + a1*st_sz, y + b10*k0)
        k2 = st_sz * func(x_n + a2*st_sz, y + b20*k0 + b21*k1)
        k3 = st_sz * func(x_n + a3*st_sz, y + b30*k0 + b31*k1 + b32*k2)
        k4 = st_sz * func(x_n + a4*st_sz, y + b40*k0 + b41*k1 + b42*k2 + b43*k3)
        k5 = st_sz * func(x_n + a5*st_sz, y + b50*k0 + b51*k1 + b52*k2 + b53*k3 + b54*k4)
        k6 = st_sz * func(x_n + a6*st_sz, y + b60*k0 + b62*k2 + b63*k3 + b64*k4 + b65*k5)
        # Getting to the slope is the whole point of this mess
        dy = c0*k0 + c2*k2 + c3*k3 + c4*k4 + c5*k5
        # Determine the estimated change in slope by comparing the output coefficients for each RK coefficient
        E = (c0 - d0)*k0 + (c2 - d2)*k2 + (c3 - d3)*k3 + (c4 - d4)*k4 + (c5 - d5)*k5 - d6*k6
        # Find the estimated error using a sum of squares method
        e = math.sqrt(np.sum(E**2)/len(y))
        # we don't know if the new value i
        hNext = 0.9*st_sz*(tol/e)**0.2
        pcnt = (i/iter_lim)*100
        psolv = (x_n/x_f)*100
        print('Correction limit : {:1.2f}% x-domain solved: {:1.2f}%'.format(pcnt, psolv))
        # If approximated error is within tolerance, accept this integration step and move on
        if e <= tol:
            # Store the new result
            i = i-1
            y = y + dy
            # Increment the x-value by the new step size
            x_n = x_n + st_sz
            # Add the new values into the output vector
            X.append(x_n)
            Y.append(y)
            # Check to break the loop when we have reached the desired x-value
            if stopper == 1: break # Reached end of x-range
            # Set limits on how much the next step size can increase to avoid missing data points
            if abs(hNext) > 10*abs(st_sz): hNext = 10*st_sz
            # Determine if the algorithm has reached the end of the dataset
            if (st_sz > 0.0) == ((x_n + hNext) >= x_f):
                hNext = x_f - x_n
                # Sets the break condition for the next loop iteration
                stopper = 1
                print('Success! Reached the end of the data set.')
            # Setting k0 to k6 * (next step size) / (current step size) forces the algorithm to use the 4th order formula for the next step
            k0 = k6*hNext/st_sz
        else:
            # The error estimate is outside the required threshold to move on, we need to redo the calculation with a smaller step size
            if abs(hNext) < abs(st_sz)*0.1 : hNext = st_sz*0.1
            # Set up k0 to go through the 5th order RK method on the next iteration because the error was no good.
            k0 = k0*hNext/st_sz
        # Set the next iteration step size
        st_sz = hNext

    # Returns the arrays for x and y values
    return np.array(X), np.array(Y)


def example_ode45py():
    '''
    Three linked bungee jumpers are depicted in figure P25.26.
    If the bungee cords are adealized as linear springs,
    the following equations based on force balances can be developed:

    m_1 * d2x/dt2 = m_1*g + k2*(x_2 - x_1) - k1*x_1
    m_2 * d2x/dt2 = m_2*g + k3*(x_3 - x_2) - k2*(x_1 - x_2)
    m_3 * d2x/dt2 = m_3*g + k3*(x_2 - x_3)

    Where:
    m_i = mass of jumper i (kg)
    kj = spring for bungee cord j (N/m)
    x_i = displacement of jumper i starting from the top of the fall and positive downward (m)
    g = gravity = 9.81 m/s**2

    Solve these equations for the positions and velocities of the three
    jumpers given the initial conditions that all positions and velocities are zero at t=0.
    Use the following parameters for your calculations:

    m_1 = 60 kg
    m_2 = 70 kg
    m_3 = 80 kg
    k1 = k3 = 50 N/m
    k2 = 100 N/m
    '''
    def fcn(t, x):
        # Instantiate an array of zeros to hold the array elements
        fcn = np.zeros(6)
        # Mass of each jumper
        m1 = 60 # kg
        m2 = 70 # kg
        m3 = 80 # kg
        # Gravity
        g = 9.81 # m/s**2
        # Bungee Cord Spring         k1 = 50 # (N/m)
        k2 = 100 # (N/m)
        k3 = 50 # (N/m)
        # Velocities for each of the jumpers:
        # dx_1/dt = v1
        fcn[0] = x[3]
        # dx_2/dt = v2
        fcn[1] = x[4]
        # dx_3/dt = v3
        fcn[2] = x[5]
        # Displacements for each jumper:
        # dv_1/dt = g + (k2/m1)*x_2 - ((k1 + k2)/m1)*x_1
        fcn[3] = g + (k2/m1)*x[1] - ((k1 + k2)/m1)*x[0]
        # dv_2/dt = g + (k2/m2)*x_1 - ((k1 + k2)/m2)*x_2 + (k3/m2)*x_3
        fcn[4] = g + (k2/m2)*x[0] - ((k2 + k3)/m2)*x[1] + (k3/m2)*x[2]
        # dv_3/dt = g + (k3/m3)*x_2 - (k3/m3)*x_3
        fcn[5] = g + (k3/m3)*x[1] - (k3/m3) * x[2]
        return fcn
    # Array for start and ending times
    t = np.array([0, 100])
    # Array holding the initial values (velocity and displacement) for each of the jumpers (everything starts at zero)
    x = np.array([0]*6)
    # Feed ode45py almost exactly like you would in MATLAB
    X, Y = ode45py(fcn, t, x, iter_lim=200000)

    # Displacement data is stored in the first three columns
    plt.plot(X, Y[:,0], label='m1 = 60kg')
    plt.plot(X, Y[:,1], label='m2 = 70kg')
    plt.plot(X, Y[:,2], label='m3 = 80kg')
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.title('Bungee Jumper Displacement')
    plt.legend()
    plt.show()

    # Velocity data is stored in the last three columns
    plt.plot(X, Y[:,3], label='m1 = 60kg')
    plt.plot(X, Y[:,4], label='m2 = 70kg')
    plt.plot(X, Y[:,5], label='m3 = 80kg')
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (m/s)')
    plt.title('Bungee Jumper Velocity')
    plt.legend()
    plt.show()
