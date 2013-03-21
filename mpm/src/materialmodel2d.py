import numpy as np

class Error(Exception):
    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg    

class JacobianError(Error):
    def __init__(self, expr, msg):
        Error.__init__(self,expr,msg)

#===============================================================================
class MaterialModel:
    # Defines material models - accessed using getStress 
    #  - actual computation done in static methods   
    # Returns stress tensor and jacobian of deformation
    def __init__(self, modelName):
        self.modelName = modelName               # Selects Material Model
        
    def getStress( self, props, F ):
        model = getattr( self, self.modelName )
        S,Ja = model(props, F);    
        return (S,Ja)
    

    @staticmethod
    def planeStrainNeoHookean( props, F ):
        # Props - poisson, E
        I2 = np.eye(2)
        v = props['poisson']
        E = props['modulus']
        l = E * v / ((1.+v)*(1.-2.*v))
        m = 0.5 * E / (1.+v)
        Ja = np.linalg.det(F)
        if Ja < 0: 
            raise JacobianError('planeStrainNeoHookean', 'Negative Jacobian')
        S = I2*l*np.log(Ja)/Ja + m/Ja * (np.dot(F, F.T) - I2)
        return (S,Ja)
    
    @staticmethod
    def planeStrainNeoHookeanMaxStress( props, F ):
        # Props - poisson, E
        I2 = np.eye(2)
        v = props['poisson']
        E = props['modulus']
        sMax = props['maxStress']
        l = E * v / ((1.+v)*(1.-2.*v))
        m = 0.5 * E / (1.+v)
        Ja = np.linalg.det(F)
        if Ja < 0: print Ja
        S = I2*l*np.log(Ja)/Ja + m/Ja * (np.dot(F, F.T) - I2)
        if vonMises(S) > sMax: S = I2*0.
        return (S,Ja)
    
    @staticmethod
    def vonMises( S ):
        return np.sqrt( S[0,0]*S[0,0] - S[0,0]*S[1,1] + S[1,1]*S[1,1] +
            3.*S[1,0]*S[0,1] )