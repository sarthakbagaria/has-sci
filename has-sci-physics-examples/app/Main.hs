module Main where

import Physics.Quantum.Algebra
import Physics.Quantum.Examples.HarmonicOscillator

import Graphics.Gloss         (Picture)


smallBoxConfig :: OscillatorInABoxConfig
smallBoxConfig = OscillatorInABoxConfig { plankConstantConfig = 6.62607004e-16 -- scaling length dimension by 10^9 so that 1 corresponds to 1 nm
                                        , particleMassConfig = 10e-31 -- roughly mass of an electron
                                        , oscillatorFrequencyConfig = 1000
                                        , halfWidthConfig = 5
                                        , halfHeightConfig = 5
                                        , initialAmplitudeConfig = gaussianWave 0 0 1
                                        , fpsConfig = 5
                                        , timeFactorConfig = 0.1
                                        }



main :: IO ()
main  = tupleToArgs simulateOscillation $! sparseModelSimulator smallBoxConfig
--main = tupleToArgs simulateOscillation $! arrayModelSimulator smallBoxConfig



tupleToArgs :: (VectorSpace k v, SimpleOperator k v o, NormalizeVector v)
            => ((TimeStep -> o) -> v  -> (v -> Picture) -> OscillatorInABoxConfig -> IO ())
            -> (( (TimeStep -> o),  v,  (v -> Picture), OscillatorInABoxConfig ) -> IO () )
tupleToArgs f (op, vec, dis, conf) = f op vec dis conf
