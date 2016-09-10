{-# LANGUAGE TypeOperators, MultiParamTypeClasses #-}

module Physics.Classical.Examples where

import Physics.Classical

import qualified Numeric.Units.Dimensional as D
import Numeric.Units.Dimensional                    ((*~), (/~))
import Numeric.Units.Dimensional.SIUnits            (newton, meter, kilo, gram, second)

import Control.Monad            (forM_)

import Graphics.Gloss
import Graphics.Gloss.Data.ViewPort                 (ViewPort)



------------
-- Particles with gravity

data PointParticle = PointParticle { pointParticleMass :: D.Mass Float
                                   , pointParticlePosition :: Vec D.DLength Float
                                   , pointParticleVelocity :: Vec (D.DLength D./ D.DTime) Float
                                   }

instance HasMass PointParticle Float where
    mass particle = pointParticleMass particle

instance HasCenterOfGravity PointParticle Float  where
    centerOfGravity particle = pointParticlePosition particle

instance HasVelocity PointParticle Float where
     velocity particle = pointParticleVelocity particle

instance ModifiableCenterOfGravity PointParticle Float where
     changeCenterOfGravity newCoG particle = particle { pointParticlePosition = newCoG }

instance ModifiableVelocity PointParticle Float where
     changeVelocity newVelocity particle = particle { pointParticleVelocity = newVelocity }



window :: Display
window = InWindow "Classical Gravity Simulator" (500, 500) (250, 250)

background :: Color
background = black

initialModel :: [ PointParticle ]
initialModel = [ PointParticle { pointParticleMass = 10e8 *~ (kilo gram)
                               , pointParticlePosition = Vec (120 *~ meter, 0 *~ meter)
                               , pointParticleVelocity = Vec (0 *~ (meter D./ second), 10 *~ (meter D./ second))
                               }
               , PointParticle { pointParticleMass = 10e12 *~ (kilo gram)
                               , pointParticlePosition = Vec (0 *~ meter, 180 *~ meter)
                               , pointParticleVelocity = Vec ((-5) *~ (meter D./ second), 0 *~ (meter D./ second))
                               }
               , PointParticle { pointParticleMass = 10e13 *~ (kilo gram)
                               , pointParticlePosition = Vec (0 *~ meter, 0 *~ meter)
                               , pointParticleVelocity = Vec (0 *~ (meter D./ second), 0 *~ (meter D./ second))
                               }
               ]
  
drawing :: [PointParticle] -> Picture
drawing particles = pictures $ map (\particle -> translate ((fst $ getVec $ pointParticlePosition particle) /~ meter)
                                                           ((snd $ getVec $ pointParticlePosition particle) /~ meter) $
                                                     color white $
                                                         circleSolid $ log ((pointParticleMass particle) /~ kilo gram)
                                   ) particles
                                                 

fps :: Int
fps = 60

simulation :: ViewPort -> Float -> [PointParticle] -> [PointParticle]
simulation _ time oldStates = nBodyGravityStep oldStates $ time *~ second

    
simulateGravity :: IO ()
simulateGravity = simulate window
                           black
                           fps
                           initialModel
                           drawing
                           simulation
