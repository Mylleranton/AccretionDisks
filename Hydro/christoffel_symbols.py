import numpy as np
#### Christoffel symbols ####
#### Schwartzschild metrics ####


##
## Radial - radial component
## Theta - angle
## a - black hole spin
## i,j,k - christoffel indices 
## 0 == Time, 1 == Radial, 2 == angular (theta), 3 == azimuthal (phi)
def christoffel(radial, theta, a, i, j, k): 
    if i == 0:
        # 0?? ################################
        if j == 0:
            #00? -------------------------------
            if k == 0:
                #000
                return 0
            elif k == 1:
                #001
                return ( 1./((-2.+radial)*radial) ) 
            elif k == 2:
                #002
                return 0
            elif k == 3:
                #003
                return 0
            else:
                print('ERROR: 00?')
                return 0
        elif j == 1:
            #01? -------------------------------
            if k == 0:
                #010
                return ( 1./((-2.+radial)*radial) )
            elif k == 1:
                #011
                return 0
            elif k == 2:
                #012
                return 0
            elif k == 3:
                #013
                return 0
            else:
                print('ERROR: 01?')
                return 0
        elif j == 2:
            #02? -------------------------------
            if k == 0:
                #020
                return 0
            elif k == 1:
                #021
                return 0
            elif k == 2:
                #022
                return 0
            elif k == 3:
                #023
                return 0
            else:
                print('ERROR: 02?')
                return 0
        elif j == 3:
            #03? -------------------------------
            if k == 0:
                #030
                return 0
            elif k == 1:
                #031
                return 0
            elif k == 2:
                #032
                return 0
            elif k == 3:
                #033
                return 0
            else:
                print('ERROR: 03?')
                return 0
        else:
            print('ERROR: 0??')
            return 0
    elif i == 1:
        # 1?? ################################
        if j == 0:
            #10? -------------------------------
            if k == 0:
                #100
                return (  (-2.+radial)/(np.power(radial, 3.))  )
            elif k == 1:
                #101
                return 0
            elif k == 2:
                #102
                return 0
            elif k == 3:
                #103
                return 0
            else:
                print('ERROR: 10?')
                return 0
        elif j == 1:
            #11? -------------------------------
            if k == 0:
                #110
                return 0
            elif k == 1:
                #111
                return (  1./(2.*radial-np.power(radial, 2.))  )
            elif k == 2:
                #112
                return 0
            elif k == 3:
                #113
                return 0
            else:
                print('ERROR: 11?')
                return 0
        elif j == 2:
            #12? -------------------------------
            if k == 0:
                #120
                return 0
            elif k == 1:
                #121
                return 0
            elif k == 2:
                #122
                return (2.-radial)
            elif k == 3:
                #123
                return 0
            else:
                print('ERROR: 12?')
                return 0
        elif j == 3:
            #13? -------------------------------
            if k == 0:
                #130
                return 0
            elif k == 1:
                #131
                return 0
            elif k == 2:
                #132
                return 0
            elif k == 3:
                #133
                return ( -((-2.+radial)*(np.power(np.sin(theta), 2.))) )
            else:
                print('ERROR: 13?')
                return 0
        else:
            print('ERROR: 1??')
            return 0       
    elif i == 2:
        # 2?? ################################
        if j == 0:
            #20? -------------------------------
            if k == 0:
                #200
                return 0
            elif k == 1:
                #201
                return 0
            elif k == 2:
                #202
                return 0
            elif k == 3:
                #203
                return 0
            else:
                print('ERROR: 20?')
                return 0
        elif j == 1:
            #21? -------------------------------
            if k == 0:
                #210
                return 0
            elif k == 1:
                #211
                return 0
            elif k == 2:
                #212
                return ( 1./radial )
            elif k == 3:
                #213
                return 0
            else:
                print('ERROR: 21?')
                return 0
        elif j == 2:
            #22? -------------------------------
            if k == 0:
                #220
                return 0
            elif k == 1:
                #221
                return (1./radial)
            elif k == 2:
                #222
                return 0
            elif k == 3:
                #223
                return 0
            else:
                print('ERROR: 22?')
                return 0
        elif j == 3:
            #23? -------------------------------
            if k == 0:
                #230
                return 0
            elif k == 1:
                #231
                return 0
            elif k == 2:
                #232
                return 0
            elif k == 3:
                #233
                return (-(np.cos(theta)*np.sin(theta)))
            else:
                print('ERROR: 23?')
                return 0
        else:
            print('ERROR: 2??')
            return 0
    elif i == 3:
        # 3?? ################################
        if j == 0:
            #30? -------------------------------
            if k == 0:
                #300
                return 0
            elif k == 1:
                #301
                return 0
            elif k == 2:
                #302
                return 0
            elif k == 3:
                #303
                return 0
            else:
                print('ERROR: 30?')
                return 0
        elif j == 1:
            #31? -------------------------------
            if k == 0:
                #310
                return 0
            elif k == 1:
                #311
                return 0
            elif k == 2:
                #312
                return 0
            elif k == 3:
                #313
                return (1./radial)
            else:
                print('ERROR: 31?')
                return 0
        elif j == 2:
            #32? -------------------------------
            if k == 0:
                #320
                return 0
            elif k == 1:
                #321
                return 0
            elif k == 2:
                #322
                return 0
            elif k == 3:
                #323
                return (1./(np.tan(theta)))
            else:
                print('ERROR: 32?')
                return 0
        elif j == 3:
            #33? -------------------------------
            if k == 0:
                #330
                return 0
            elif k == 1:
                #331
                return (1./radial)
            elif k == 2:
                #332
                return (1./(np.tan(theta)))
            elif k == 3:
                #333
                return 0
            else:
                print('ERROR: 33?')
                return 0
        else:
            print('ERROR: 3??')
            return 0
    else:
        print('ERROR: This should not happen... If you see this, go home for the day')
        return 0
