import numpy as np
import matplotlib.pyplot as plt 
import datetime, os

g = 9.81 # Gravity
rho = 1.2 # Density of air at sea level
# coeffDrag = 0.5 # Drag Coefficient

initialVelocity = 860 # m/s
initialAngle = np.radians(30) # degrees
barrelTwistRate = 12 # Inches/rotation
initialSpin = (12/barrelTwistRate)*(initialVelocity*(2*np.pi*3.28084))


# M80 FMJ Ball
r = 0.00381 # m Bullet radius (7.62x51 mm NATO) 0.3085in
ogiveLength = 0.0194818 # m 0.767in
ogiveRadius = 0.0762 # m 3.00in
meplatDiameter = 0.001524 # m 0.060in
boatTailLength = 0.004445 # m .175in
boatTailAngle = np.radians(9.0) # deg
oal = 0.028956 # m 1.140in
bulletMass = 9.52544 # kg 9.52544g = 147gr

# Assume homogenous density and approximate mass distribution with cones and cylinders
diameter = 2*r
areaSection = np.pi*r**2
BaseRadius = r - boatTailLength/np.tan(boatTailAngle)
volumeBaseConeSmall = np.pi*BaseRadius**3*np.tan(boatTailAngle)/3
volumeBaseConeBig = np.pi*r**3*np.tan(boatTailAngle)/3
volumeOgive = np.pi*r**2*ogiveLength/3
volumeBoatTail = (np.pi*np.tan(boatTailAngle)/3)*(r**3 - BaseRadius**3)
volumeBearing = np.pi*r**2*(oal - boatTailLength - ogiveLength)
volumeTotal = volumeOgive + volumeBoatTail + volumeBearing
density = bulletMass/1000*volumeTotal # <--------------------------- Modified for Expanding Bullet Mass
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



def get_filename(prefix, suffix, base_path):
    '''
    Gets a unique file name in the base path.
    
    Appends date and time information to file name and adds a number
    if the file name is stil not unique.

    prefix = Homework assignment name

    suffix = Extension

    base_path = Location of log file
    '''
    # Set base filename for compare
    fileNameBase = base_path + prefix + "_" + datetime.datetime.now().strftime("%y%m%d_%H%M%S")
    # Set base for numbering system if filename exists
    num = 1
    # Generate complete filename to check existence
    fileName = fileNameBase + suffix
    # Find a unique filename
    while os.path.isfile(fileName):
        # if the filename is not unique, add a number to the end of it
        fileName = fileNameBase + "_" + str(num) + suffix
        # increments the number in case the filename is still not unique
        num = num + 1
    return fileName


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
    # Cash-Karp coefficients for RK algorithm -
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
        e = np.sqrt(np.sum(E**2)/len(y))
        # we don't know if the new value i
        hNext = 0.9*st_sz*(tol/e)**0.2
        pcnt = (i/iter_lim)*100
        psolv = (x_n/x_f)*100
        print('Correction limit : {:1.2f}  x-domain solved: {:1.2f}#'.format(pcnt, psolv))
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


def bulletSim():

    def fcn(t, x):
        # Instantiate an array of zeros to hold the array elements
        fcn = np.zeros(4)

        fcn[0] = Ux = x[2]
        fcn[1] = Uy = x[3]
        Cd = 0.6366198/(density*r^2*(x[2]^2+x[1]^2))
        fcn[2] = -Cd*Ux/bulletMass
        fcn[3] = -g - Cd*Uy/bulletMass

        return fcn
    # Array for start and ending times
    t = np.array([0, 60])
    x = np.array([0]*4)
    x[0] = 0
    x[1] = 0
    x[2] = initialVelocity*np.cos(initialAngle)
    x[3] = initialVelocity*np.sin(initialAngle)

    # Feed ode45py almost exactly like you would in MATLAB
    X, Y = ode45py(fcn, t, x, tol=1.0e-8, iter_lim=20000)
    dataFile = get_filename('trivialData','.csv','./log_files/')
    heading = 'time(s), x_pos, y_pos, x_vel, y_vel'

    np.savetxt(dataFile, (X, Y[:,0], Y[:,1], Y[:,2], Y[:,3]), delimiter=',', header=heading, newline='\n')

    # plotFile1 = get_filename('Trajectory','.png','./log_files/')
    # # plt.plot(X, Y[:,0], label='X Position')
    # # plt.plot(X, Y[:,1], label='Y Position')
    # # plt.plot(X, Y[:,2], label='X Velocity')
    # # plt.plot(X, Y[:,3], label='Y Velocity')
    # plt.plot(Y[:,0],Y[:,1], label='Displacement')
    # plt.xlabel('Displacement X (m)')
    # plt.ylabel('Displacement Y (m)')
    # plt.title('Trajectory Data')
    # plt.legend()
    # # plt.savefig(plotFile1, bbox_inches='tight')
    # plt.show()

    plotFile2 = get_filename('velocityAndTime','.png','./log_files/')

    vel = np.sqrt(Y[:,2]**2 + Y[:,3]**2)
    for k, v in enumerate(vel):
        if v < 343:
            soundIndex = k
            soundValue = v
            break
    plt.plot(X, vel)
    plt.scatter(X[soundIndex], soundValue, label='Speed of Sound' ,color='r')
    plt.annotate('({:1.2f} s)'.format(X[soundIndex]), (X[soundIndex] + 2 , soundValue))
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (m/s)')
    plt.title('Velocity Data')
    plt.legend()
    # plt.savefig(plotFile2, bbox_inches='tight')
    plt.show()

    plotFile3 = get_filename('Trajectory_withSound','.png','./log_files/')
    plt.plot(Y[:,0],Y[:,1])
    plt.scatter(Y[soundIndex,0], Y[soundIndex,1], label='Sound Barrier' ,color='r')
    plt.annotate('({:1.0f}, {:1.0f})'.format(Y[soundIndex,0], Y[soundIndex,1]), (Y[soundIndex,0]-400, Y[soundIndex,1]+0), ha='right')
    plt.xlabel('Displacement X (m)')
    plt.ylabel('Displacement Y (m)')
    plt.title('Trajectory Data')
    plt.legend()
    # plt.savefig(plotFile3, bbox_inches='tight')
    plt.show()
bulletSim()