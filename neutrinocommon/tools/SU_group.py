import numpy as np

i = np.complex(0,1)

def Conmutator(X,Y):
    return np.dot(X,Y) - np.dot(Y,X)

def AntiConmitator(X,Y):
    return np.dot(X,Y) + np.dot(Y,X)

def GellMannMatrices():
    
    lambda1 = [[0,1,0],[1,0,0],[0,0,0]]
    lambda2 = [[0,-i,0],[i,0,0],[0,0,0]]
    lambda3 = [[1,0,0],[0,-1,0],[0,0,0]]
    lambda4 = [[0,0,1],[0,0,0],[1,0,0]]
    lambda5 = [[0,0,-i],[0,0,0],[i,0,0]]
    lambda6 = [[0,0,0],[0,0,1],[0,1,0]]
    lambda7 = [[0,0,0],[0,0,-i],[0,i,0]]
    lambda8 = [[1/np.sqrt(3),0,0],[0,0,0],[0,0,-2/np.sqrt(3)]]
    
    return [lambda1,lambda2,lambda3,lambda4,
            lambda5,lambda8,lambda7,lambda8]
    
if __name__ == '__main__':
    
    lam =  GellMannMatrices()
    
    print np.dot(lam[1],lam[1])
    
    print Conmutator(lam[0],lam[1])
    
    