module Main where

import Physics.Quantum.Algebra
import Physics.Quantum.Examples.HarmonicOscillator

import Graphics.Gloss         (Picture)


smallBoxConfig :: OscillatorInABoxConfig
smallBoxConfig = OscillatorInABoxConfig { plankConstantConfig = 6.62607004e-16 -- scaling length dimension by 10^9 so that 1 corresponds to 1 nm
                                        , particleMassConfig = 10e-31 -- roughly mass of an electron
                                        , oscillatorFrequencyConfig = 10000
                                        , halfWidthConfig = 10
                                        , halfHeightConfig = 10
                                        , initialAmplitudeConfig = gaussianWave 0 0 2
                                        , fpsConfig = 30
                                        , timeFactorConfig = 0.5
                                        }

bigBoxConfig :: OscillatorInABoxConfig
bigBoxConfig = OscillatorInABoxConfig { plankConstantConfig = 6.62607004e-16 -- scaling length dimension by 10^9 so that 1 corresponds to 1 nm
                                      , particleMassConfig = 10e-25
                                      , oscillatorFrequencyConfig = 10e7
                                      , halfWidthConfig = 100
                                      , halfHeightConfig = 100
                                      , initialAmplitudeConfig = gaussianWave 30 30 30
                                      , fpsConfig = 1
                                      , timeFactorConfig = 0.01
                                      }
  


main :: IO ()
main  = tupleToArgs simulateOscillation $! sparseMatModelSimulator smallBoxConfig
--main = tupleToArgs simulateOscillation $! arrayModelSimulator smallBoxConfig



tupleToArgs :: (VectorSpace k v, SimpleOperator k v o, NormalizeVector v)
            => ((TimeStep -> o) -> v  -> (v -> Picture) -> OscillatorInABoxConfig -> IO ())
            -> (( (TimeStep -> o),  v,  (v -> Picture), OscillatorInABoxConfig ) -> IO () )
tupleToArgs f (op, vec, dis, conf) = f op vec dis conf
