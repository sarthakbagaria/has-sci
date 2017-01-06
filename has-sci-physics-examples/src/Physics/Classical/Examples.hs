{-# LANGUAGE TypeOperators, MultiParamTypeClasses, FlexibleContexts #-}

module Physics.Classical.Examples where

import Physics.Classical

import qualified Numeric.Units.Dimensional as D
import qualified Numeric.Units.Dimensional.Quantities as D
import Numeric.Units.Dimensional                    ((*~), (/~))
import Numeric.Units.Dimensional.SIUnits            (newton, meter, kilo, gram, second, coulomb, nano)

import Control.Monad            (forM_)

import Graphics.Gloss
import Graphics.Gloss.Data.ViewPort                 (ViewPort)
import Graphics.Gloss.Raster.Field                  (rgb', makePicture)



------------
-- Particles with gravity

data PointParticle = PointParticle { pointParticleMass :: D.Mass Float
                                   , pointParticlePosition :: Vec D.DLength Float
                                   , pointParticleVelocity :: Vec (D.DLength D./ D.DTime) Float
                                   , pointParticleElectricCharge :: D.Quantity D.DElectricCharge Float
                                   }

instance HasMass PointParticle Float where
    mass particle = pointParticleMass particle

instance HasElectricCharge PointParticle Float where
    electricCharge particle = pointParticleElectricCharge particle

instance HasCenterOfGravity PointParticle Float  where
    centerOfGravity particle = pointParticlePosition particle

instance HasVelocity PointParticle Float where
     velocity particle = pointParticleVelocity particle

instance ModifiableCenterOfGravity PointParticle Float where
     changeCenterOfGravity newCoG particle = particle { pointParticlePosition = newCoG }

instance ModifiableVelocity PointParticle Float where
     changeVelocity newVelocity particle = particle { pointParticleVelocity = newVelocity }



    
simulateGravity :: IO ()
simulateGravity = simulate window
                           black
                           fps
                           initialModel
                           drawing
                           simulation
    where
        window :: Display
        window = InWindow "Classical Gravity Simulator" (500, 500) (250, 250)

        background :: Color
        background = black

        initialModel :: [ PointParticle ]
        initialModel = [ PointParticle { pointParticleMass = 10e8 *~ (kilo gram)
                                       , pointParticlePosition = Vec (120 *~ meter, 0 *~ meter)
                                       , pointParticleVelocity = Vec (0 *~ (meter D./ second), 10 *~ (meter D./ second))
                                       , pointParticleElectricCharge = 0 *~ coulomb 
                                       }
                       , PointParticle { pointParticleMass = 10e12 *~ (kilo gram)
                                       , pointParticlePosition = Vec (0 *~ meter, 180 *~ meter)
                                       , pointParticleVelocity = Vec ((-5) *~ (meter D./ second), 0 *~ (meter D./ second))
                                       , pointParticleElectricCharge = 0 *~ coulomb 
                                       }
                       , PointParticle { pointParticleMass = 10e13 *~ (kilo gram)
                                       , pointParticlePosition = Vec (0 *~ meter, 0 *~ meter)
                                       , pointParticleVelocity = Vec (0 *~ (meter D./ second), 0 *~ (meter D./ second))
                                       , pointParticleElectricCharge = 0 *~ coulomb
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




simulateElectricCharges :: IO ()
simulateElectricCharges = simulate window
                                   black
                                   fps
                                   initialModel
                                   drawing
                                   simulation
    where
        window :: Display
        window = InWindow "Classical Electric Simulator" (500, 500) (250, 250)

        background :: Color
        background = black

        initialModel :: [ PointParticle ]
        initialModel = [
                        -- three small particles with small charges
                        PointParticle { pointParticleMass = 1e5 *~ (kilo gram)
                                      , pointParticlePosition = Vec (100 *~ nano meter, (-60) *~ nano meter)
                                      , pointParticleVelocity = Vec (0 *~ (nano meter D./ second), 10 *~ (nano meter D./ second))
                                      , pointParticleElectricCharge = 1e-16 *~ coulomb 
                                      }
                       , PointParticle { pointParticleMass = 1e5 *~ (kilo gram)
                                       , pointParticlePosition = Vec ((-150) *~ nano meter, 20 *~ nano meter)
                                       , pointParticleVelocity = Vec (10 *~ (nano meter D./ second), 0 *~ (nano meter D./ second))
                                       , pointParticleElectricCharge = 1e-16 *~ coulomb
                                       }
                       , PointParticle { pointParticleMass = 1e5 *~ (kilo gram)
                                       , pointParticlePosition = Vec (150 *~ nano meter, 80 *~ nano meter)
                                       , pointParticleVelocity = Vec ((-5) *~ (nano meter D./ second), (-5) *~ (nano meter D./ second))
                                       , pointParticleElectricCharge = 1e-16 *~ coulomb
                                       }

                       -- four massive particles with high charges, dominating the electric potential and remaining stationary
                       , PointParticle { pointParticleMass = 1e20 *~ (kilo gram)
                                       , pointParticlePosition = Vec (250 *~ (nano meter), 250 *~ (nano meter))
                                       , pointParticleVelocity = Vec (0 *~ (meter D./ second), 0 *~ (meter D./ second))
                                       , pointParticleElectricCharge = 1e-10 *~ coulomb
                                       }
                       , PointParticle { pointParticleMass = 1e20 *~ (kilo gram)
                                       , pointParticlePosition = Vec ((-250) *~ nano meter, 250 *~ nano meter)
                                       , pointParticleVelocity = Vec (0 *~ (meter D./ second), 0 *~ (meter D./ second))
                                       , pointParticleElectricCharge = 1e-10 *~ coulomb
                                       }
                       , PointParticle { pointParticleMass = 1e20 *~ (kilo gram)
                                       , pointParticlePosition = Vec ((-250) *~ (nano meter), (-250) *~ (nano meter))
                                       , pointParticleVelocity = Vec (0 *~ (meter D./ second), 0 *~ (meter D./ second))
                                       , pointParticleElectricCharge = 1e-10 *~ coulomb
                                       }
                       , PointParticle { pointParticleMass = 1e20 *~ (kilo gram)
                                       , pointParticlePosition = Vec (250 *~ (nano meter), (-250) *~ (nano meter))
                                       , pointParticleVelocity = Vec (0 *~ (meter D./ second), 0 *~ (meter D./ second))
                                       , pointParticleElectricCharge = 1e-10 *~ coulomb
                                       }
                       ]

        drawing :: [PointParticle] -> Picture
        drawing particles = pictures $ potentialChart :
                                       ( map (\particle -> translate ((fst $ getVec $ pointParticlePosition particle) /~ nano meter)
                                                                     ((snd $ getVec $ pointParticlePosition particle) /~ nano meter) $
                                                               color white $
                                                                   circleSolid $ (log $ ((pointParticleMass particle) /~ kilo gram) + 10) / 5
                                             ) particles
                                       )


        -- We are computing the background potential only once for the entire simulation
        -- This is because the dominant charges are very massive and will not move much in the simulation
        -- So the potential will not change much
        -- To compute the potential dynamically for each frame, put this computation inside drawing function and use current particle state
        -- instead of initial model
        -- Also, we are computing max potential and min potential approximately based on heuristics, and not exactly
        -- This is to tune the contours and colors by hand,
        -- and also for dynamic potential computing min and max exactly will add quite a lot to computation
        potentialChart = makePicture 500 500 1 1 ( \(x,y) -> let potential = electricPotentialNormalized (x*250, y*250)
                                                                 -- 50 countour lines spaced at frequency 1/10
                                                                 (red, green, blue) = case (round $! 500 * potential) `mod` 10 == 0 of
                                                                     True
                                                                         |                      potential < -1    -> (0,0,0)
                                                                         |potential >= -1    && potential < -0.75 -> (0, potential*4 + 4, 0)
                                                                         |potential >= -0.75 && potential < 0     -> (0, 0, potential*4/3 + 1)
                                                                         |potential >= 0     && potential < 1     -> (potential, 0, 0)
                                                                         |potential >= 1                          -> (1,0,0)
                                                                     False -> (0,0,0)
                                                     
                                                             in rgb' red green blue
                                                 )
            where
                electricPotentialSI (x,y) = (nBodyElectricPotential initialModel $ Vec (x *~ nano meter, y *~ nano meter)) /~ (newton D.* meter D./ coulomb)
                electricPotentialNormalized pos = {-(electricPotentialSI pos) / maxPotential-} (electricPotentialSI pos - potentialRangeCenter) / ((maxPotential - minPotential)/2)
            
                -- This computes the potential due to largest charge about 100 nano meters away from it
                -- This potential is high enough for our simulation to not consider variation of potential above it
                maxPotential = maximum $ map (\particle -> abs $ ( electricFieldConstant D.* (electricCharge particle) D./ (100 *~ nano meter) ) /~ (newton D.* meter D./ coulomb)
                                             ) initialModel
                -- a rough approximation to lower bound for potential for the given simulation
                minPotential = ( electricFieldConstant D.* (foldr (D.+) D._0 $ map electricCharge initialModel) D./ (500 *~ nano meter) ) /~ (newton D.* meter D./ coulomb)/2
                potentialRangeCenter = (minPotential + maxPotential)


        fps :: Int
        fps = 60

        simulation :: ViewPort -> Float -> [PointParticle] -> [PointParticle]
        simulation _ time oldStates = nBodyElectricStep oldStates $ time *~ second
