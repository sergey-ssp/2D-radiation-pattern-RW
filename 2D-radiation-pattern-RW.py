import numpy as np
import matplotlib.pyplot as plt

def radiation_pattern(a, b, L, lamda, SMPL_x, SMPL_y, dx, dy):
    # calculate wavenumber "k" :
    k = 2. * np.pi / lamda

    # calculate overall dimensions of the scan plane :
    scene_x = (SMPL_x - 1) * dx
    scene_y = (SMPL_y - 1) * dy

    # initiate matrices for the distributions of electromagnetic field at the scan plane :
    G = [ [0 for ix in range(SMPL_x)] for jy in range(SMPL_y) ]
    Gsq = [ [0 for ix in range(SMPL_x)] for jy in range(SMPL_y) ]

    # initiate matrices for the distributions of electromagnetic field at the middle axes of the scan plane :
    G_midX = [0 for jy in range(SMPL_y)]
    G_midY = [0 for ix in range(SMPL_x)]

    # calculate the EM field distribution at scan plane (Fraunhofer region) :
    for ix in range(SMPL_x) :
        xx = -scene_x / 2 + ix * dx

        for jy in range(SMPL_y) :
            yy = -scene_y / 2 + jy * dy
            var1 = 4.*np.pi*a*(L**3) / k / yy
            var1 = var1 / ( (np.pi*L)**2 + (a*k*xx)**2 )
            var2 = np.cos(0.5*k*a*xx/L)
            var3 = np.sin(0.5*k*b*yy/L)
            G[jy][ix] = var1 * var2 * var3
            Gsq[jy][ix] = G[jy][ix] **2

            # store the distribution at the middle axes :
            if ix == SMPL_x / 2 :
                G_midX[jy] = Gsq[jy][ix]

            if jy == SMPL_y / 2 :
                G_midY[ix] = Gsq[jy][ix]

    # plot output distributions :
    x = np.linspace(-scene_x/2, scene_x/2, SMPL_x, endpoint=True)
    y = np.linspace(-scene_y/2, scene_y/2, SMPL_y, endpoint=True)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure()

    fig.suptitle("Intensity of EM Field $G^2$ at L = {}".format(L), fontsize=14)

    ax1 = fig.add_subplot(121)
    extent = [0 , SMPL_x, 0 , SMPL_y]
    pattern = ax1.imshow(Gsq, cmap='jet', extent=[-scene_x/2, scene_x/2, -scene_y/2, scene_y/2])
    fig.colorbar(pattern)
    ax1.set_title('Pattern')
    ax1.set_xlabel('X (mm)')
    ax1.set_ylabel('Y (mm)')

    ax2 = fig.add_subplot(222)
    ax2.plot(x, G_midY, label='y = 0')
    ax2.legend()
    ax2.set_xlabel('X (mm)')
    ax2.set_ylabel('$G^2$')

    ax3 = fig.add_subplot(224)
    ax3.plot(x, G_midX, label='x = 0')
    ax3.legend()
    ax3.set_xlabel('Y (mm)')
    ax3.set_ylabel('$G^2$')

    plt.show()




if __name__ == "__main__":
    # dimensions of rectangular waveguide 'a' x 'b' :
    a = 3.6  #sinusoidal distribution (x)
    b = 1.8  #uniform distribution (y)

    # distance from waveguide output to scan plane 'L' :
    L = 50.

    # radiation wavelenght 'lamda' :
    lamda = 2.15

    # dimensions of output matrix :
    SMPL_x = 500
    SMPL_y = 500

    # scales 'dx' and 'dy' :
    dx = 0.2
    dy = 0.2

    # calculation is valid only for D^2 / lamda / L << 1
    radiation_pattern(a, b, L, lamda, SMPL_x, SMPL_y, dx, dy)
