import numpy as np

class MaterialModel:
    # Defines material models - accessed using getStress 
    #  - actual computation done in static methods    
    def __init__(self, modelName):
        self.modelName = modelName               # Selects Material Model
        
    def getStress( self, props, F ):
        model = getattr( self, self.modelName )
        S,Ja = model(props, F);    
        return (S,Ja)

    @staticmethod
    def planeStrainNeoHookean( props, F ):
        I2 = np.eye(2)
        v = props.poisson
        l = props.E * v / ((1.+v)*(1.-2.*v))
        m = 0.5 * props.E / (1.+v)
        Ja = np.linalg.det(F)
        S = I2*l*np.log(Ja)/Ja + m/Ja * ( np.dot( F, F.T ) - I2 )
        return (S,Ja)