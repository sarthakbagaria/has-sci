module Main where

import Physics.Quantum.Algebra
import Physics.Quantum.Examples.HarmonicOscillator
import Physics.Classical.Examples                   (simulateGravity, simulateElectricCharges)

import Graphics.Gloss         (Picture)
import System.Environment     (getArgs)


-- quantum oscillation simulator
main :: IO ()
main = do
    args <- getArgs
    case args of
        "classical_gravity" : _ -> simulateGravity
        "classical_electric_potential" : _ -> simulateElectricCharges
        "quantum_oscillation" : _ -> tupleToArgs simulateOscillation $! sparseMatModelSimulator smallBoxConfig
                                     -- tupleToArgs simulateOscillation $! sparseMatModelSimulator bigBoxConfig
                                     -- tupleToArgs simulateOscillation $! sparseModelSimulator smallBoxConfig
        -- simulate electric charges by default
        _ -> simulateElectricCharges




tupleToArgs :: (VectorSpace k v, SimpleOperator k v o, NormalizeVector v)
            => ((TimeStep -> o) -> v  -> (v -> Picture) -> OscillatorInABoxConfig -> IO ())
            -> (( (TimeStep -> o),  v,  (v -> Picture), OscillatorInABoxConfig ) -> IO () )
tupleToArgs f (op, vec, dis, conf) = f op vec dis conf
