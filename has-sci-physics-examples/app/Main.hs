module Main where

import Physics.Quantum.Algebra
import Physics.Quantum.Examples.HarmonicOscillator
import Physics.Classical.Examples

import Graphics.Gloss         (Picture)


smallBoxConfig :: OscillatorInABoxConfig
smallBoxConfig = OscillatorInABoxConfig { plankConstantConfig = 6.62607004e-16 -- scaling length dimension by 10^9 so that 1 corresponds to 1 nm
                                        , particleMassConfig = 10e-31 -- roughly mass of an electron
                                        , oscillatorFrequencyConfig = 10000
                                        , halfWidthConfig = 10
                                        , halfHeightConfig = 10
                                          -- tune variance of wave based on mass and frequency (and Plank's constant) to get a nice simulation
                                        , initialAmplitudeConfig = gaussianWave 0 0 3.8845e-6
                                        , fpsConfig = 30
                                        , timeFactorConfig = 0.5
                                        }

bigBoxConfig :: OscillatorInABoxConfig
bigBoxConfig = OscillatorInABoxConfig { plankConstantConfig = 6.62607004e-16 -- scaling length dimension by 10^9 so that 1 corresponds to 1 nm
                                      , particleMassConfig = 10e-31
                                      , oscillatorFrequencyConfig = 10000
                                      , halfWidthConfig = 100
                                      , halfHeightConfig = 100
                                        -- tune variance of wave based on mass and frequency (and Plank's constant) to get a nice simulation
                                      , initialAmplitudeConfig = gaussianWave 0 0 3.8845e-6
                                      , fpsConfig = 20
                                      , timeFactorConfig = 0.5
                                      }  

  


main :: IO ()
main = tupleToArgs simulateOscillation $! sparseMatModelSimulator smallBoxConfig
--main = tupleToArgs simulateOscillation $! sparseMatModelSimulator bigBoxConfig
--main = simulateGravity
--main = tupleToArgs simulateOscillation $! sparseModelSimulator smallBoxConfig
--main = tupleToArgs simulateOscillation $! arrayModelSimulator smallBoxConfig



tupleToArgs :: (VectorSpace k v, SimpleOperator k v o, NormalizeVector v)
            => ((TimeStep -> o) -> v  -> (v -> Picture) -> OscillatorInABoxConfig -> IO ())
            -> (( (TimeStep -> o),  v,  (v -> Picture), OscillatorInABoxConfig ) -> IO () )
tupleToArgs f (op, vec, dis, conf) = f op vec dis conf
