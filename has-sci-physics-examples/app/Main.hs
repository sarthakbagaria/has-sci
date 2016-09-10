module Main where

import Physics.Quantum.Algebra
import Physics.Quantum.Examples.HarmonicOscillator
import Physics.Classical.Examples

import Graphics.Gloss         (Picture)



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
