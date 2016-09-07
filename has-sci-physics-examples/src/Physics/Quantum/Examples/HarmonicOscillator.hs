{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}

module Physics.Quantum.Examples.HarmonicOscillator where

import Physics.Quantum.Algebra

import Numeric.Additive.Class          (Additive((+)), Abelian(..))
import Numeric.Algebra.Class           ( Multiplicative((*), pow1p)
                                       , Monoidal(..), Semiring(..), LeftModule(..), RightModule(..)
                                       )

import Numeric.Algebra.Division        (Division(..))
import Numeric.Domain.Euclidean        (Euclidean(..))
import Numeric.Domain.PID              (PID(..))
import Numeric.Domain.UFD              (UFD(..))
import Numeric.Domain.GCD              (GCDDomain(..))
import Numeric.Domain.Integral         (IntegralDomain(..))
import Numeric.Algebra.Unital          (Unital(..))
import Numeric.Algebra.Unital.UnitNormalForm (UnitNormalForm(..))
import Numeric.Decidable.Zero          (DecidableZero(..))
import Numeric.Semiring.ZeroProduct    (ZeroProductSemiring(..))
import Numeric.Algebra.Commutative     (Commutative(..))
import Numeric.Decidable.Units         (DecidableUnits(..))
import Numeric.Algebra.Unital.UnitNormalForm  (UnitNormalForm)
import Numeric.Ring.Class              (Ring(..))
import Numeric.Decidable.Associates    (DecidableAssociates(..))
import Numeric.Rig.Class               (Rig(..))
import Numeric.Additive.Group          (Group(..))

import GHC.Natural                     (Natural)
import GHC.Float                       (double2Float)

import Data.Array.Repa.Eval            (Elt(..))
import Data.Vector.Unboxed             (Unbox(..))

  
import Data.Complex                    (Complex(..), magnitude)

import Prelude hiding ((+),(-),(/),(*), (^))
import qualified Prelude as P 

import qualified Data.Map as M
import qualified Data.Array.Repa as R
import Data.Array.Repa                              (Z(..), (:.)(..))
import Data.Functor.Identity (runIdentity)

import Graphics.Gloss
import Graphics.Gloss.Data.ViewPort                 (ViewPort)

import Data.Traversable      (for)
import Data.Foldable         (foldl')

import Debug.Trace           (trace)
import Control.Concurrent    (threadDelay)



--------------------
-- Add dimensions to quantities
-- Since we are writing operators as endomorphisms of vector spaces
-- if we have to give dimensiones to operators/eigenvalues,
-- we may need to declare the algebra of dimensions as instance of field too
-- Or we may need to relax the condition of operator to homomorphims of vector spaces
-------------------

--------------------
-- This simulator is super-inefficient at the moment
-- only runs on a lattice of around 3x3
-- Use textSimulateOscillation
-- as 3x3 pixels on graphical simulation may not be very visible 
-------------------


-------------------
-- For field instance of Complex Double

instance Euclidean (Complex Double) where degree _ = Just 1
instance Division (Complex Double)
instance PID (Complex Double)
instance UFD (Complex Double)
instance GCDDomain (Complex Double)
instance IntegralDomain (Complex Double)
instance Unital (Complex Double) where one = 1 :+ 0
instance DecidableZero (Complex Double) where isZero = ((0 :+ 0) ==)
instance ZeroProductSemiring (Complex Double)
instance Commutative (Complex Double)
instance Multiplicative (Complex Double) where (*) = (P.*)
instance DecidableUnits (Complex Double) where recipUnit a = case a of
                                                (1 :+ 0) -> Just a
                                                (-1 :+ 0) -> Just a
                                                _ -> Nothing
instance UnitNormalForm (Complex Double)
instance Ring (Complex Double)
instance Monoidal (Complex Double) where zero = 0 :+ 0
instance Semiring (Complex Double)
instance Additive (Complex Double) where (+) = (P.+)
instance DecidableAssociates (Complex Double) where isAssociate = (==)
instance Rig (Complex Double)
instance Group (Complex Double)
instance LeftModule Natural (Complex Double) where
    n .* m = (fromIntegral $ toInteger n) P.* m
instance RightModule Natural (Complex Double) where
    m *. n = m P.* (fromIntegral n)
instance LeftModule Integer (Complex Double) where
    n .* m = (fromIntegral $ toInteger n) P.* m
instance RightModule Integer (Complex Double) where
    m *. n = m P.* (fromIntegral n)
instance Abelian (Complex Double)


-- For use by repa
instance Elt (Complex Double) where
    touch (a :+ b) = touch a >> touch b 



--------------

type Position = (Int, Int)

positionRange :: [Position]
positionRange = [(x,y) | x <- [(-halfWidth) .. halfWidth], y <- [-(halfHeight) .. halfHeight]]


-------------
-- Using index functions

{-
posXOp :: InfOperator (Complex Double) Position
posXOp = InfOperator $ \(x,y) -> M.fromList [((x,y), (fromIntegral x) :+ 0) ]
 
posYOp :: InfOperator (Complex Double) Position
posYOp = InfOperator $ \(x,y) -> M.fromList [((x,y), (fromIntegral y) :+ 0)]

momXOp :: InfOperator (Complex Double) Position
-- taking dx as 1
momXOp = InfOperator $ \(x,y) -> M.fromList [ ( (x P.- 1, y), (0 :+ (plankConstant)) ),
                                              ( (x P.+ 1, y), (0 :+ (-(plankConstant))) )
                                            ]

momYOp :: InfOperator (Complex Double) Position
-- taking dy as 1
momYOp = InfOperator $ \(x,y) -> M.fromList [ ( (x,y P.- 1), 0 :+ (plankConstant) ),
                                              ( (x,y P.+ 1), 0 :+ (-(plankConstant)) )
                                            ]

scOp :: (Complex Double) -> InfOperator (Complex Double) Position
scOp scalar = InfOperator $! \pos -> M.fromList [ (pos, scalar) ]


oscillationOperator :: Double -> InfOperator (Complex Double) Position
oscillationOperator dt = trace "OperationEnding" $! (scOp $! 1 :+ 0) + ((scOp $! 0 :+ (0 P.- (dt P./ plankConstant))) *
                                                   (
                                                     ( (momXOp `pow1p` 2) + (momYOp `pow1p` 2) * (scOp $ (1 P./ (2 P.* particleMass)) :+ 0) ) +
                                                     ( (posXOp `pow1p` 2) + (posYOp `pow1p` 2) * (scOp $ ((particleMass P./ 2) P.* oscillatorFrequency P.^ 2) :+ 0) )
                                                   )
                                                )

initialModel :: InfKet (Complex Double) Position
initialModel = InfKet $! \pos -> if pos == initialParticlePosition then 1 else 0

drawing :: InfKet (Complex Double) Position -> Picture
drawing (InfKet coeffMap) = trace "Drawing" $
  pictures $ map (\(x,y) -> translate (fromIntegral x) (fromIntegral y) $!
                            color (withAlpha (double2Float $! normSquared $! coeffMap (x,y)) red) $!
                            circleSolid 1
                 ) $ pixels
    where pixels = positionRange
          normSquared cn = (magnitude cn) P.^ 2

-}



------------------
-- Using arrays and matrices


-- Helper functions to translate between two dimensional position and one dimensional array

indexRange :: Int
indexRange = (2 * halfWidth + 1) * (2 * halfHeight + 1)

positionToIndex :: Position -> Int
positionToIndex pos = (fst pos + halfWidth) * (2 * halfHeight + 1) + (snd pos + halfHeight)

indexToPosition :: Int -> Position
indexToPosition index = ( (floor ((fromIntegral index) P./ (fromIntegral $! 2 * halfHeight + 1))) - halfWidth
                        , (index `mod` (2 * halfHeight + 1)) - halfHeight
                        )


-- position and momentum operators

posXOp :: LatticeOperator (Complex Double) Position
posXOp = LatticeOperator $! R.fromFunction (Z :. indexRange :. indexRange) (\(Z:.i:.j) -> if i==j
                                                                                          then (fromIntegral $! fst $! indexToPosition i) :+ 0
                                                                                          else 0
                                                                           )
  
posYOp :: LatticeOperator (Complex Double) Position
posYOp = LatticeOperator $! R.fromFunction (Z :. indexRange :. indexRange) (\(Z:.i:.j) -> if i==j
                                                                                          then (fromIntegral $! snd $! indexToPosition i) :+ 0
                                                                                          else 0
                                                                           )
  
momXOp :: LatticeOperator (Complex Double) Position
momXOp = LatticeOperator $! R.fromFunction (Z :. indexRange :. indexRange) (\(Z:.i:.j) -> let xo = fst $! indexToPosition j
                                                                                              yo = snd $! indexToPosition j
                                                                                              xn = fst $! indexToPosition i
                                                                                              yn = snd $! indexToPosition i
                                                                                          in
                                                                                            case yn==yo of
                                                                                              False -> 0
                                                                                              True
                                                                                                -- taking dx as 1
                                                                                                |xn == xo P.+ 1 -> 0 :+ (plankConstant)
                                                                                                |xn == xo P.- 1 -> 0 :+ (-(plankConstant))
                                                                                                |otherwise -> 0
                                                                           )

momYOp :: LatticeOperator (Complex Double) Position
momYOp = LatticeOperator $! R.fromFunction (Z :. indexRange :. indexRange) (\(Z:.i:.j) -> let xo = fst $! indexToPosition j
                                                                                              yo = snd $! indexToPosition j
                                                                                              xn = fst $! indexToPosition i
                                                                                              yn = snd $! indexToPosition i
                                                                                          in
                                                                                            case xn==xo of
                                                                                              False -> 0
                                                                                              True
                                                                                                -- taking dy as 1
                                                                                                |yn == yo P.+ 1 -> 0 :+ (plankConstant)
                                                                                                |yn == yo P.- 1 -> 0 :+ (-(plankConstant))
                                                                                                |otherwise -> 0
                                                                           )

-- a scalar operator or scalar multiplication
-- can make operators instances of Module and remove this
scOp :: (Complex Double) -> LatticeOperator (Complex Double) Position
scOp scalar = LatticeOperator $! R.fromFunction (Z :. indexRange :. indexRange) (\(Z:.i:.j) -> case i==j of
                                                                                                   False -> 0
                                                                                                   True -> scalar
                                                                                )
  

oscillationOperator :: Double -> LatticeOperator (Complex Double) Position
oscillationOperator dt = trace "OperationEnding" $!
      (scOp $! 1 :+ 0) + ( (scOp $! 0 :+ (0 P.- (dt P./ plankConstant))) *
                           (
                             ( ( (momXOp `pow1p` 2) + (momYOp `pow1p` 2) ) * (scOp $ (1 P./ (2 P.* particleMass)) :+ 0) ) +
                             ( ( (posXOp `pow1p` 2) + (posYOp `pow1p` 2) ) * (scOp $ ((particleMass P./ 2) P.* oscillatorFrequency P.^ 2) :+ 0) )
                           )
                         )

initialModel :: LatticeKet (Complex Double) Position
initialModel = LatticeKet $! R.fromFunction (Z :. indexRange) (\(Z:.i) -> if ( i == positionToIndex initialParticlePosition )
                                                                          then 1
                                                                          else 0
                                                              )


drawing :: LatticeKet (Complex Double) Position -> Picture
drawing (LatticeKet vector) = trace "Drawing" $
  pictures $ map (\(x,y) -> translate (fromIntegral x) (fromIntegral y) $!
                            color (withAlpha (double2Float $ normSquared $! (vector R.! ( Z:.(positionToIndex (x,y))))) red) $!
                            circleSolid 1
                 ) $ pixels
    where pixels = positionRange
          normSquared cn = (magnitude cn) P.^ 2



textSimulateOscillation :: IO ()
textSimulateOscillation = do
    return $! trace "Start Calculating Operator." ()
    operator <- return $ oscillationOperator $! (1 P./ (fromIntegral fps)) P.* timeFactor
    print $! runIdentity $! R.computeUnboxedP $! unLatticeOperator $! operator
    return $! trace "Finish Calculating Operator." ()
    go operator initialModel

    where
      go :: LatticeOperator (Complex Double) Position -> LatticeKet (Complex Double) Position -> IO ()
      go operator ket@(LatticeKet vector) = do
              _ <- return $! runIdentity $! R.computeUnboxedP $! vector
              print "Done"
              threadDelay 1000000
              -- need to normalize the vector
              let nextVector = operator `op` ket 
              go operator nextVector


-----------------------------------------
-- final simulator function
-- common to funcion and array methods

simulateOscillation :: IO ()
simulateOscillation = do
    -- computing the operator
    return $! trace "Start Calculating Operator." ()
    operator <- return $ oscillationOperator $! (1 P./ (fromIntegral fps)) P.* timeFactor
    return $! trace "Finish Calculating Operator." ()
    simulate window
             black
             fps
             initialModel
             drawing
             (simulation operator)
    where
        -- need to normalize the vector after operator action
        simulation operator _ _ ket = trace "Simulating" $! operator `op` ket



------------------------
-- Configuration

-- scaling length dimension by 10^9 so that 1 corresponds to i nm
plankConstant :: Double
plankConstant = 6.62607004e-16

-- roughly mass of an electron
particleMass :: Double
particleMass = 10e-31

oscillatorFrequency :: Double
oscillatorFrequency = 100

halfWidth :: Int
halfWidth = 1

halfHeight :: Int
halfHeight = 1

window :: Display
window = InWindow "Nice Window" (2*halfWidth, 2*halfHeight) (halfWidth, halfHeight)

background :: Color
background = black

fps :: Int
fps = 1

initialParticlePosition :: Position
initialParticlePosition = (1,0)

-- simulation will play 10 times slower than real time
-- mainly to keep time steps small
-- while keeping fps less
timeFactor :: Double
timeFactor = 0.01
