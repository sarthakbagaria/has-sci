{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}

module Physics.Quantum.Examples.HarmonicOscillator where

import Physics.Quantum.Algebra

import Numeric.Additive.Class          (Additive((+)), Abelian)
import Numeric.Algebra.Class           ( Multiplicative((*), pow1p)
                                       , Monoidal(..), Semiring, LeftModule(..), RightModule(..)
                                       )

import Numeric.Algebra.Division        (Division(..))
import Numeric.Domain.Euclidean        (Euclidean(..))
import Numeric.Domain.PID              (PID(..))
import Numeric.Domain.UFD              (UFD)
import Numeric.Domain.GCD              (GCDDomain(..))
import Numeric.Domain.Integral         (IntegralDomain(..))
import Numeric.Algebra.Unital          (Unital(..))
import Numeric.Algebra.Unital.UnitNormalForm (UnitNormalForm(..))
import Numeric.Decidable.Zero          (DecidableZero(..))
import Numeric.Semiring.ZeroProduct    (ZeroProductSemiring)
import Numeric.Algebra.Commutative     (Commutative)
import Numeric.Decidable.Units         (DecidableUnits(..))
import Numeric.Ring.Class              (Ring(..))
import Numeric.Decidable.Associates    (DecidableAssociates(..))
import Numeric.Rig.Class               (Rig(..))
import Numeric.Additive.Group          (Group(..))

import GHC.Natural                     (Natural)
import GHC.Float                       (double2Float)

import Data.Array.Repa.Eval            (Elt(touch))

  
import Data.Complex                    (Complex(..), magnitude, realPart)

import Prelude hiding ((+),(-),(/),(*), (^))
import qualified Prelude as P 

import qualified Data.Map as M
import qualified Data.Array.Repa as R
import Data.Array.Repa                              (Z(..), (:.)(..))
import Data.Functor.Identity (runIdentity)

import qualified Data.Eigen.SparseMatrix as ES
import qualified Data.Eigen.Matrix as ESD
import Data.Eigen.Matrix                            (CComplex)
import Foreign.C.Types                              (CDouble)

import Data.Maybe                                   (catMaybes)
import qualified Data.List as L

import Graphics.Gloss
import Graphics.Gloss.Data.ViewPort                 (ViewPort)
import Graphics.Gloss.Raster.Field                  (rgb', makePicture)

import Debug.Trace           (trace)
import Control.Concurrent    (threadDelay)



--------------------
-- Add dimensions to quantities
-- Since we are writing operators as endomorphisms of vector spaces
-- if we have to give dimensiones to operators/eigenvalues,
-- we may need to declare the algebra of dimensions as instance of field too
-- Or we may need to relax the condition of operator to homomorphims of vector spaces
-------------------




------------------------------------------
-- Configure
------------------------------------------

data OscillatorInABoxConfig = OscillatorInABoxConfig { plankConstantConfig :: Double
                                                     , particleMassConfig :: Double
                                                     , oscillatorFrequencyConfig :: Double
                                                     , halfWidthConfig :: Int
                                                     , halfHeightConfig :: Int
                                                     , initialAmplitudeConfig :: Position -> Complex Double
                                                     , fpsConfig :: Int
                                                     -- the simulation will play this times faster than real time
                                                     -- mainly to keep time steps small while keeping fps low
                                                     -- the time step in each simulation will be this divided by fps
                                                     , timeFactorConfig :: TimeStep
                                                     }



type TimeStep = Double

data Position = Position { posX :: Int , posY :: Int }
    deriving (Eq, Ord, Show)


type Variance = Double


gaussianWave :: Int -> Int -> Variance -> (Position -> Complex Double)
gaussianWave xCenter yCenter var = \(Position x y) -> 
    ( exp ( (
              - ( ((fromIntegral (x - xCenter))P.^2) P.+ ((fromIntegral (y - yCenter) )P.^2) )
            ) P./ (2 P.* (var P.^ 2))
          )
    ) :+ 0



smallBoxConfig :: OscillatorInABoxConfig
smallBoxConfig = OscillatorInABoxConfig { plankConstantConfig = 6.62607004e-16 -- scaling length dimension by 10^9 so that 1 corresponds to 1 nm
                                        , particleMassConfig = 10e-31 -- roughly mass of an electron
                                        , oscillatorFrequencyConfig = 10000
                                        , halfWidthConfig = 10
                                        , halfHeightConfig = 10
                                          -- tune variance of wave based on mass and frequency (and Plank's constant) to get a nice simulation
                                        , initialAmplitudeConfig = gaussianWave 0 0 3.8845e-6
                                        , fpsConfig = 30
                                        , timeFactorConfig = 1
                                        }

bigBoxConfig :: OscillatorInABoxConfig
bigBoxConfig = OscillatorInABoxConfig { plankConstantConfig = 6.62607004e-16 -- scaling length dimension by 10^9 so that 1 corresponds to 1 nm
                                      , particleMassConfig = 10e-31
                                      , oscillatorFrequencyConfig = 10000
                                      , halfWidthConfig = 100
                                      , halfHeightConfig = 100
                                        -- tune variance of wave based on mass and frequency (and Plank's constant) to get a nice simulation
                                      , initialAmplitudeConfig = gaussianWave 0 0 3.8845e-6
                                      , fpsConfig = 30
                                      , timeFactorConfig = 1
                                      }  


---------------------
-- Helper functions to generate lattice and translate between two dimensional position and one dimensional array
-- may not need to modify these
--------------------

positionRangeConfig :: OscillatorInABoxConfig -> [Position]
positionRangeConfig config = [Position x y | x <- [-halfWidth .. halfWidth], y <- [-halfHeight .. halfHeight]]
    where halfWidth = halfWidthConfig config
          halfHeight = halfHeightConfig config


indexRangeConfig :: OscillatorInABoxConfig -> Int
indexRangeConfig config = (2 * halfWidth + 1) * (2 * halfHeight + 1)
    where halfWidth = halfWidthConfig config
          halfHeight = halfHeightConfig config


positionToIndexConfig :: OscillatorInABoxConfig -> Position -> Int
positionToIndexConfig config pos = ((posX pos + halfWidth) * (2 * halfHeight + 1)) + (posY pos + halfHeight)
    where halfWidth = halfWidthConfig config
          halfHeight = halfHeightConfig config
                              

indexToPositionConfig :: OscillatorInABoxConfig -> Int -> Position
indexToPositionConfig config index = Position ( (floor ((fromIntegral index) P./ (fromIntegral $! 2 * halfHeight + 1))) - halfWidth )
                                              ( (index `mod` (2 * halfHeight + 1)) - halfHeight )
    where halfWidth = halfWidthConfig config
          halfHeight = halfHeightConfig config
  




------------------------------------
-- simulate oscillation using Eigen library's sparse matrices
------------------------------------


sparseMatModelSimulator :: OscillatorInABoxConfig
                        -> ( TimeStep -> SparseMatOperator (Complex Double) (CComplex CDouble) Position
                           , SparseMatKet (Complex Double) (CComplex CDouble) Position
                           , SparseMatKet (Complex Double) (CComplex CDouble) Position -> Picture
                           , OscillatorInABoxConfig
                           ) 
sparseMatModelSimulator config = (oscillationOperator, initialModel, drawing, config)
    where
        -- shorthands for config params
        positionRange = positionRangeConfig config
        halfWidth = halfWidthConfig config
        halfHeight = halfHeightConfig config
        indexRange = indexRangeConfig config

        positionToIndex = positionToIndexConfig config
        indexToPosition = indexToPositionConfig config

        plankConstant = plankConstantConfig config
        particleMass = particleMassConfig config
        oscillatorFrequency = oscillatorFrequencyConfig config
        
        initialAmplitude = initialAmplitudeConfig config
        
      
        -- position and momentum operators
        
        posXOp :: SparseMatOperator (Complex Double) (CComplex CDouble) Position
        posXOp = SparseMatOperator $! ES.fromList indexRange indexRange $! map (\pos@(Position x _) -> (positionToIndex pos, positionToIndex pos, fromIntegral x)) positionRange

        posYOp :: SparseMatOperator (Complex Double) (CComplex CDouble) Position
        posYOp = SparseMatOperator $! ES.fromList indexRange indexRange $! map (\pos@(Position _ y) -> (positionToIndex pos, positionToIndex pos, fromIntegral y)) positionRange

        momXOp :: SparseMatOperator (Complex Double) (CComplex CDouble) Position
        -- taking dx as 2
        -- <x'| P_x |a>  =  -ih ( d/dx(a)|_x' ) = - ( a(x'+1) - a(x'-1) ) * ih/2   = - (<x'+1|a> - <x'-1|a>) * ih / 2
        momXOp = SparseMatOperator $! ES.fromList indexRange indexRange $! ( (map (\pos@(Position x y) -> ( positionToIndex pos
                                                                                                          , positionToIndex (Position (x P.+ 1) y)
                                                                                                          , 0 :+ (- plankConstant/2)
                                                                                                          )
                                                                                  )
                                                                              -- removing max x from position range to keep operator within lattice
                                                                              -- that would give 0 contribution anyway
                                                                              [Position i j | i <- [(-halfWidth) .. (halfWidth-1)], j <- [-(halfHeight) .. halfHeight]]
                                                                             )
                                                                             ++
                                                                             (map (\pos@(Position x y) -> ( positionToIndex pos
                                                                                                          , positionToIndex (Position (x P.- 1) y)
                                                                                                          , 0 :+ (  plankConstant/2)
                                                                                                          )
                                                                                  )
                                                                              -- removing min x from position range to keep operator within lattice
                                                                              -- that would give 0 contribution anyway
                                                                              [Position i j | i <- [-(halfWidth-1) .. (halfWidth)], j <- [-(halfHeight) .. halfHeight]]
                                                                             )
                                                                           )
        
        momYOp :: SparseMatOperator (Complex Double) (CComplex CDouble) Position
        -- taking dy as 2
        -- <y'| P_y |a>  =  -ih ( d/dy(a)|_y' ) = - ( a(y'+1) - a(y'-1) ) * ih/2   = - (<y'+1|a> - <y'-1|a>) * ih / 2
        momYOp = SparseMatOperator $! ES.fromList indexRange indexRange $! ( (map (\pos@(Position x y) -> ( positionToIndex pos
                                                                                                          , positionToIndex (Position x (y P.+ 1))
                                                                                                          , 0 :+ (- plankConstant/2)
                                                                                                          )
                                                                                  )
                                                                              -- removing max y from position range to keep operator within lattice
                                                                              -- that would give 0 contribution anyway
                                                                              [Position i j | i <- [(-halfWidth) .. (halfWidth)], j <- [-(halfHeight) .. (halfHeight-1)]]
                                                                             )
                                                                             ++
                                                                             (map (\pos@(Position x y) -> ( positionToIndex pos
                                                                                                          , positionToIndex (Position x (y P.- 1))
                                                                                                          , 0 :+ (  plankConstant/2)
                                                                                                          )
                                                                                  )
                                                                              -- removing min y from position range to keep operator within lattice
                                                                              -- that would give 0 contribution anyway
                                                                              [Position i j | i <- [(-halfWidth) .. (halfWidth)], j <- [-(halfHeight-1) .. halfHeight]]
                                                                             )
                                                                           )

        momSquaredOp = (momXOp `pow1p` 2) + (momYOp `pow1p` 2)
        posSquaredOp = (posXOp `pow1p` 2) + (posYOp `pow1p` 2)
          

        idOp :: SparseMatOperator (Complex Double) (CComplex CDouble) Position
        idOp = SparseMatOperator $! ES.fromList indexRange indexRange $! map (\pos@(Position _ _) -> (positionToIndex pos, positionToIndex pos, 1)) positionRange


        oscillationOperator :: Double -> SparseMatOperator (Complex Double) (CComplex CDouble) Position
        oscillationOperator dt =
              idOp + ( (0 :+ (0 P.- (dt P./ plankConstant))) .*
                       (
                         ( (((particleMass P./ 2) P.* (oscillatorFrequency P.^ 2)) :+ 0) .* posSquaredOp ) +
                         ( ((1 P./ (2 P.* particleMass)) :+ 0)                           .* momSquaredOp )
                       )
                     )

          
        initialModel :: SparseMatKet (Complex Double) (CComplex CDouble) Position
        initialModel = normalizeVector $! SparseMatKet $! ES.fromList indexRange 1 $! map ( \pos -> ( positionToIndex pos
                                                                                                    , 0
                                                                                                    , initialAmplitude pos
                                                                                                    )
                                                                                          ) positionRange

        -- bitmap drawing
        drawing :: SparseMatKet (Complex Double) (CComplex CDouble) Position -> Picture
        drawing (SparseMatKet mat) =
            makePicture (2*halfWidth) (2*halfHeight)
                        1 1
                        ( \(x,y) -> let index = positionToIndex $! Position (round $! x P.* (P.fromIntegral halfWidth))
                                                                            (round $! y P.* (P.fromIntegral halfHeight))
                                        amp = mat ES.! (index, 0)
                                        den = double2Float $! normSquared amp
                                    in {-trace (show $ (index, amp, den)) $ -} rgb' den 0 0
                        )
            where
                normSquared cn = (magnitude cn) P.^ 2


        {-
        -- drawing with objects
        drawing :: SparseMatKet (Complex Double) (CComplex CDouble) Position -> Picture
        drawing (SparseMatKet mat) = {- trace "Drawing" $ -}
            pictures $ catMaybes $ map (\elem -> case elem of
                                           (index, 0, amp) -> let pos = indexToPosition index
                                                              in Just $! translate (fromIntegral $! posX pos) (fromIntegral $! posY pos) $!
                                                                 color (withAlpha (double2Float $! normSquared amp) red) $!
                                                                 circleSolid 1
                                           _ -> trace "OutOfBounds: " Nothing
                                       ) $! ES.toList $! mat
          where
                normSquared cn = (magnitude cn) P.^ 2
        -}





------------------------------------
-- simulate oscillation using sparse maps
------------------------------------


sparseModelSimulator :: OscillatorInABoxConfig
                     -> ( TimeStep -> SparseOperator (Complex Double) Position
                        , SparseKet (Complex Double) Position
                        , SparseKet (Complex Double) Position -> Picture
                        , OscillatorInABoxConfig
                        ) 
sparseModelSimulator config = (oscillationOperator, initialModel, drawing, config)
    where
        -- shorthands for config params
        positionRange = positionRangeConfig config
        halfWidth = halfWidthConfig config
        halfHeight = halfHeightConfig config

        plankConstant = plankConstantConfig config
        particleMass = particleMassConfig config
        oscillatorFrequency = oscillatorFrequencyConfig config
        
        initialAmplitude = initialAmplitudeConfig config
        
      
        -- position and momentum operators
        
        posXOp :: SparseOperator (Complex Double) Position
        posXOp = SparseOperator $! M.fromList $! map (\pos@(Position x _) -> ((pos,pos), ((fromIntegral x) :+ 0))) positionRange

        posYOp :: SparseOperator (Complex Double) Position
        posYOp = SparseOperator $! M.fromList $! map (\pos@(Position _ y) -> ((pos,pos), ((fromIntegral y) :+ 0))) positionRange

        momXOp :: SparseOperator (Complex Double) Position
        -- taking dx as 2
        -- <x'| P_x |a>  =  -ih ( d/dx(a)|_x' ) = - ( a(x'+1) - a(x'-1) ) * ih/2   = - (<x'+1|a> - <x'-1|a>) * ih / 2        
        momXOp = SparseOperator $! M.fromList $! ( (map (\pos@(Position x y) -> ( (pos, Position (x P.+ 1) y), 0 :+ (- plankConstant/2) ))
                                                        -- removing max x from position range to keeo operator within lattice
                                                        -- that would give 0 contribution anyway
                                                        [Position i j | i <- [(-halfWidth) .. (halfWidth-1)], j <- [-(halfHeight) .. halfHeight]]
                                                   )
                                                   ++
                                                   (map (\pos@(Position x y) -> ( (pos, Position (x P.- 1) y), 0 :+ (  plankConstant/2) ))
                                                        -- removing min x from position range to keeo operator within lattice
                                                        -- that would give 0 contribution anyway
                                                        [Position i j | i <- [-(halfWidth-1) .. (halfWidth)], j <- [-(halfHeight) .. halfHeight]]
                                                   )
                                                 )

        momYOp :: SparseOperator (Complex Double) Position
        -- taking dy as 2
        -- <y'| P_y |a>  =  -ih ( d/dy(a)|_y' ) = - ( a(y'+1) - a(y'-1) ) * ih/2   = - (<y'+1|a> - <y'-1|a>) * ih / 2
        momYOp = SparseOperator $! M.fromList $! ( (map (\pos@(Position x y) -> ( (pos, Position x (y P.+ 1)), 0 :+ (- plankConstant/2) ))
                                                        -- removing max y from position range to keeo operator within lattice
                                                        -- that would give 0 contribution anyway
                                                        [Position i j | i <- [(-halfWidth) .. (halfWidth)], j <- [-(halfHeight) .. (halfHeight-1)]]
                                                   )
                                                   ++
                                                   (map (\pos@(Position x y) -> ( (pos, Position x (y P.- 1)), 0 :+ (  plankConstant/2) ))
                                                        -- removing min y from position range to keeo operator within lattice
                                                        -- that would give 0 contribution anyway
                                                        [Position i j | i <- [(-halfWidth) .. (halfWidth)], j <- [-(halfHeight-1) .. halfHeight]]
                                                   )
                                                                              
                                                 )

        momSquaredOp = (momXOp `pow1p` 2) + (momYOp `pow1p` 2)
        posSquaredOp = (posXOp `pow1p` 2) + (posYOp `pow1p` 2)
          

        idOp :: SparseOperator (Complex Double) Position
        idOp = SparseOperator $! M.fromList $! map (\pos@(Position _ _) -> ((pos,pos), 1)) positionRange


        oscillationOperator :: Double -> SparseOperator (Complex Double) Position
        oscillationOperator dt =
              idOp + ( (0 :+ (0 P.- (dt P./ plankConstant))) .*
                       (
                         ( (((particleMass P./ 2) P.* (oscillatorFrequency P.^ 2)) :+ 0) .* posSquaredOp ) +
                         ( ((1 P./ (2 P.* particleMass)) :+ 0)                           .* momSquaredOp )
                       )
                     )

          
        initialModel :: SparseKet (Complex Double) Position
        initialModel = normalizeVector $! SparseKet $! M.fromList $! map ( \pos -> (pos, initialAmplitude pos) ) positionRange

        -- bitmap drawing
        drawing :: SparseKet (Complex Double) Position -> Picture
        drawing (SparseKet coeffMap) =
            makePicture (2*halfWidth) (2*halfHeight)
                        1 1
                        ( \(x,y) -> let pos = Position (round $! x P.* (P.fromIntegral halfWidth))
                                                       (round $! y P.* (P.fromIntegral halfHeight))
                                        maybeAmp = M.lookup pos coeffMap
                                    in case maybeAmp of
                                           Nothing -> rgb' 0 0 0
                                           Just amp -> let den = double2Float $! normSquared amp
                                                       in {-trace (show $ (index, amp, den)) $ -} rgb' den 0 0
                        )
            where
                normSquared cn = (magnitude cn) P.^ 2

        {-
        -- drawing with objects
        drawing :: SparseKet (Complex Double) Position -> Picture
        drawing (SparseKet coeffMap) = {- trace "Drawing" $ -} 
            pictures $ map (\((Position x y), amp) -> translate (fromIntegral x) (fromIntegral y) $!
                                                              color (withAlpha (double2Float $! normSquared amp) red) $!
                                                                  circleSolid 1
                           ) $! M.toList $! coeffMap
            where
                normSquared cn = (magnitude cn) P.^ 2
        -}



------------------------------------
-- Using repa arrays and matrices
-- Matrix multiplication turns out to be the bottleneck
-- Probably because it's not designed for sparse matrices
-- can work for around halfWidth=2, halfHeight=2
------------------------------------


arrayModelSimulator :: OscillatorInABoxConfig
                    -> ( TimeStep -> LatticeOperator (Complex Double) Position
                       , LatticeKet (Complex Double) Position
                       , LatticeKet (Complex Double) Position -> Picture
                       , OscillatorInABoxConfig
                       )
arrayModelSimulator config = (oscillationOperator, initialModel, drawing, config)
    where
        -- shorthands for config params
        positionRange = positionRangeConfig config
        halfWidth = halfWidthConfig config
        halfHeight = halfHeightConfig config

        indexRange = indexRangeConfig config
        positionToIndex = positionToIndexConfig config
        indexToPosition = indexToPositionConfig config

        plankConstant = plankConstantConfig config
        particleMass = particleMassConfig config
        oscillatorFrequency = oscillatorFrequencyConfig config
        
        initialAmplitude = initialAmplitudeConfig config
        
          
        -- position and momentum operators
      
        posXOp :: LatticeOperator (Complex Double) Position
        posXOp = LatticeOperator $! R.fromFunction (Z :. indexRange :. indexRange) (\(Z:.i:.j) -> if i==j
                                                                                                  then (fromIntegral $! posX $! indexToPosition i) :+ 0
                                                                                                  else 0
                                                                                   )

        posYOp :: LatticeOperator (Complex Double) Position
        posYOp = LatticeOperator $! R.fromFunction (Z :. indexRange :. indexRange) (\(Z:.i:.j) -> if i==j
                                                                                                  then (fromIntegral $! posY $! indexToPosition i) :+ 0
                                                                                                  else 0
                                                                                   )

        momXOp :: LatticeOperator (Complex Double) Position
        -- taking dx as 2
        -- <x'| P_x |a>  =  -ih ( d/dx(a)|_x' ) = - ( a(x'+1) - a(x'-1) ) * ih/2   = - (<x'+1|a> - <x'-1|a>) * ih / 2        
        momXOp = LatticeOperator $! R.fromFunction (Z :. indexRange :. indexRange) (\(Z:.i:.j) -> let xo = posX $! indexToPosition j
                                                                                                      yo = posY $! indexToPosition j
                                                                                                      xn = posX $! indexToPosition i
                                                                                                      yn = posY $! indexToPosition i
                                                                                                  in
                                                                                                    case yn==yo of
                                                                                                      False -> 0
                                                                                                      True
                                                                                                        |xn == xo P.+ 1 -> 0 :+ (  plankConstant/2)
                                                                                                        |xn == xo P.- 1 -> 0 :+ (- plankConstant/2)
                                                                                                        |otherwise -> 0
                                                                                   )

        momYOp :: LatticeOperator (Complex Double) Position
        -- taking dy as 2
        -- <y'| P_y |a>  =  -ih ( d/dy(a)|_y' ) = - ( a(y'+1) - a(y'-1) ) * ih/2   = - (<y'+1|a> - <y'-1|a>) * ih / 2        
        momYOp = LatticeOperator $! R.fromFunction (Z :. indexRange :. indexRange) (\(Z:.i:.j) -> let xo = posX $! indexToPosition j
                                                                                                      yo = posY $! indexToPosition j
                                                                                                      xn = posX $! indexToPosition i
                                                                                                      yn = posY $! indexToPosition i
                                                                                                  in
                                                                                                    case xn==xo of
                                                                                                      False -> 0
                                                                                                      True
                                                                                                        |yn == yo P.+ 1 -> 0 :+ (  plankConstant/2)
                                                                                                        |yn == yo P.- 1 -> 0 :+ (- plankConstant/2)
                                                                                                        |otherwise -> 0
                                                                                   )

        momSquaredOp = (momXOp `pow1p` 2) + (momYOp `pow1p` 2)
        posSquaredOp = (posXOp `pow1p` 2) + (posYOp `pow1p` 2)


        -- identity operator
        -- can make operator an instance of Unital, but don't know any one of writing the identity matrix in repa independent of dimensions
        idOp :: LatticeOperator (Complex Double) Position
        idOp = LatticeOperator $! R.fromFunction (Z :. indexRange :. indexRange) (\(Z:.i:.j) -> case i==j of
                                                                                     False -> 0
                                                                                     True -> 1
                                                                                 )

        -- simulation operator
        oscillationOperator :: Double -> LatticeOperator (Complex Double) Position
        oscillationOperator dt =
              idOp + ( (0 :+ (0 P.- (dt P./ plankConstant))) .*
                       (
                         ( (((particleMass P./ 2) P.* (oscillatorFrequency P.^ 2)) :+ 0) .* posSquaredOp ) +
                         ( ((1 P./ (2 P.* particleMass)) :+ 0)                           .* momSquaredOp )
                       )
                     )

        -- initial ket
        initialModel :: LatticeKet (Complex Double) Position
        initialModel = LatticeKet $! R.fromFunction (Z :. indexRange) $! \(Z:.i) -> initialAmplitude $! indexToPosition i

        -- bitmap drawing
        drawing :: LatticeKet (Complex Double) Position -> Picture
        drawing (LatticeKet vector) =
            makePicture (2*halfWidth) (2*halfHeight)
                        1 1
                        ( \(x,y) -> let pos = Position (round $! x P.* (P.fromIntegral halfWidth))
                                                       (round $! y P.* (P.fromIntegral halfHeight))
                                        amp = vector R.! ( Z:.(positionToIndex pos) )
                                        den = double2Float $! normSquared amp
                                    in {-trace (show $ (pos, amp, den)) $ -} rgb' den 0 0
                        )
            where
                normSquared cn = (magnitude cn) P.^ 2


                                                                                    
        {-
        -- drawing with objects
        drawing :: LatticeKet (Complex Double) Position -> Picture
        drawing (LatticeKet vector) = {- trace "Drawing" $ -}
          pictures $ map (\(Position x y) -> translate (fromIntegral x) (fromIntegral y) $!
                                    color (withAlpha (double2Float $ normSquared $! (vector R.! ( Z:.(positionToIndex (Position x y))))) red) $!
                                    circleSolid 1
                         ) $ positionRange
            where normSquared cn = (magnitude cn) P.^ 2
        -}
     



------------------------------------
-- simulate oscillation using index functions
-- doesn't work at the moment
-- no NormalizeVector instance for InfKet
------------------------------------

indexFunctionsSimulator :: OscillatorInABoxConfig
                        -> ( TimeStep -> InfOperator (Complex Double) Position
                           , InfKet (Complex Double) Position
                           , InfKet (Complex Double) Position -> Picture
                           , OscillatorInABoxConfig
                           )
indexFunctionsSimulator config = (oscillationOperator, initialModel, drawing, config)
    where
        -- shorthands for config params
        positionRange = positionRangeConfig config
        halfWidth = halfWidthConfig config
        halfHeight = halfHeightConfig config

        plankConstant = plankConstantConfig config
        particleMass = particleMassConfig config
        oscillatorFrequency = oscillatorFrequencyConfig config
        
        initialAmplitude = initialAmplitudeConfig config

      
        -- position and momentum operators
        posXOp :: InfOperator (Complex Double) Position
        posXOp = InfOperator $ \pos@(Position x _) -> M.fromList [(pos, (fromIntegral x) :+ 0) ]

        posYOp :: InfOperator (Complex Double) Position
        posYOp = InfOperator $ \pos@(Position _ y) -> M.fromList [(pos, (fromIntegral y) :+ 0)]

        momXOp :: InfOperator (Complex Double) Position
        -- taking dx as 2
        -- <x'| P_x |a>  =  -ih ( d/dx(a)|_x' ) = - ( a(x'+1) - a(x'-1) ) * ih/2   = - (<x'+1|a> - <x'-1|a>) * ih / 2        
        momXOp = InfOperator $ \(Position x y) -> M.fromList [ ( Position (x P.- 1) y, (0 :+ (  plankConstant/2)) ),
                                                               ( Position (x P.+ 1) y, (0 :+ (- plankConstant/2)) )
                                                             ]

        momYOp :: InfOperator (Complex Double) Position
        -- taking dy as 2
        -- <y'| P_y |a>  =  -ih ( d/dy(a)|_y' ) = - ( a(y'+1) - a(y'-1) ) * ih/2   = - (<y'+1|a> - <y'-1|a>) * ih / 2        
        momYOp = InfOperator $ \(Position x y) -> M.fromList [ ( Position x (y P.- 1), 0 :+ (  plankConstant/2) ),
                                                               ( Position x (y P.+ 1), 0 :+ (- plankConstant/2) )
                                                             ]

        momSquaredOp = (momXOp `pow1p` 2) + (momYOp `pow1p` 2)
        posSquaredOp = (posXOp `pow1p` 2) + (posYOp `pow1p` 2)

        scOp :: (Complex Double) -> InfOperator (Complex Double) Position
        scOp scalar = InfOperator $! \pos -> M.fromList [ (pos, scalar) ]


        oscillationOperator :: Double -> InfOperator (Complex Double) Position
        oscillationOperator dt = (scOp $! 1 :+ 0) + ((scOp $! 0 :+ (0 P.- (dt P./ plankConstant))) *
                                                     (
                                                       ( momSquaredOp * (scOp $ (1 P./ (2 P.* particleMass)) :+ 0) ) +
                                                       ( posSquaredOp * (scOp $ ((particleMass P./ 2) P.* oscillatorFrequency P.^ 2) :+ 0) )
                                                     )
                                                    )

        initialModel :: InfKet (Complex Double) Position
        initialModel = InfKet $! initialAmplitude

        -- bitmap drawing
        drawing :: InfKet (Complex Double) Position -> Picture
        drawing (InfKet coeffMap) = {- trace "Drawing" $ -}
            makePicture (2*halfWidth) (2*halfHeight)
                        1 1
                        ( \(x,y) -> let pos = Position (round $! x P.* (P.fromIntegral halfWidth))
                                                       (round $! y P.* (P.fromIntegral halfHeight))
                                        amp = coeffMap pos
                                        den = double2Float $! normSquared amp
                                    in {-trace (show $ (pos, amp, den)) $ -} rgb' den 0 0
                        )
            where
                normSquared cn = (magnitude cn) P.^ 2

        {-
        -- drawing with objects        
        drawing :: InfKet (Complex Double) Position -> Picture
        drawing (InfKet coeffMap) = {- trace "Drawing" $ -}
            pictures $ map (\(Position x y) -> translate (fromIntegral x) (fromIntegral y) $!
                                               color (withAlpha (double2Float $! normSquared $! coeffMap (Position x y)) red) $!
                                               circleSolid 1
                           ) $ positionRange
            where normSquared cn = (magnitude cn) P.^ 2
        -}




------------------------------------------
-- Simulate
------------------------------------------


simulateOscillation :: (VectorSpace k v, SimpleOperator k v o, NormalizeVector v)
                    => (TimeStep -> o)
                    -> v
                    -> (v -> Picture)
                    -> OscillatorInABoxConfig
                    -> IO ()
simulateOscillation oscillationOperator initialModel drawing config = do
    -- return $! trace "Start Calculating Operator." ()
    -- keeping fps constant to avoid recalculation of operator
    operator <- return $! oscillationOperator $! timeFactor P./ (fromIntegral fps)
    -- return $! trace "Finish Calculating Operator." ()
    simulate (InWindow "Quantum Oscillation Simulator" (2*halfWidth, 2*halfHeight) (halfWidth, halfHeight))
             black
             fps
             initialModel
             drawing
             (simulation operator)

    where
        fps = fpsConfig config
        halfHeight = halfHeightConfig config
        halfWidth = halfWidthConfig config
        timeFactor = timeFactorConfig config
        
        simulation :: (VectorSpace k v, SimpleOperator k v o, NormalizeVector v)
                   => o -> ViewPort -> Float
                   -> v -> v
        simulation operator _ _ ket = {- trace "Simulating" $! -} normalizeVector $! operator `op` ket


                        
textSimulateOscillation :: (VectorSpace k v, Show v, SimpleOperator k v o, NormalizeVector v)
                        => (TimeStep -> o)
                        -> v
                        -> OscillatorInABoxConfig
                        -> IO ()
textSimulateOscillation oscillationOperator initialModel config = do
    -- return $! trace "Start Calculating Operator." ()
    operator <- return $ oscillationOperator $! timeFactor P./ (fromIntegral fps)
    -- return $! trace "Finish Calculating Operator." ()
    go operator initialModel

    where
        fps = fpsConfig config
        timeFactor = timeFactorConfig config
        
        go :: ( VectorSpace k v, SimpleOperator k v o
              , Show v, NormalizeVector v
              )
            => o -> v -> IO ()
        go operator ket = do
                print ket
                threadDelay 1000000
                let nextMap = operator `op` ket 
                go operator $! normalizeVector nextMap





--------------------------------
--Class Instances
--------------------------------


class NormalizeVector a where
    normalizeVector :: a -> a


instance NormalizeVector (SparseMatKet (Complex Double) (CComplex CDouble) Position) where
    normalizeVector ket@(SparseMatKet mat) = (1 / (norm :+ 0)) *^ ket
        where
            norm :: Double
            -- norm = realPart $! ES.norm mat
            -- normalizing with max of norms instead of sum of squares of norm to make the states more visible in simulation
            norm = ESD.maxCoeff $ ESD.convert normSquared $ ES.toMatrix mat
            
            normSquared cn = (magnitude cn) P.^ 2


instance NormalizeVector (SparseKet (Complex Double) Position) where
    normalizeVector ket@(SparseKet ketMap) = (1 / (norm :+ 0)) *^ ket
        where
            -- ketNormSquared = M.foldl' (P.+) 0 $! M.map (\coeff -> normSquared coeff) ketMap
            -- normalizing with max of norms instead of sum of squares of norm to make the states more visible in simulation
            maxNormSquared = M.foldl' (\a b -> max a b) 0 $! M.map (\coeff -> normSquared coeff) ketMap
            normSquared cn = (magnitude cn) P.^ 2
            norm :: Double
            norm = sqrt $! maxNormSquared

        
instance NormalizeVector (LatticeKet (Complex Double) Position) where
    normalizeVector ket@(LatticeKet array) = (1 / (norm :+ 0)) *^ ket
        where
            -- ketNormSquared = R.foldlAllS (P.+) 0 $! R.map (\coeff -> normSquared coeff) array
            -- normalizing with max of norms instead of sum of squares of norm to make the states more visible in simulation
            maxNormSquared = R.foldAllS (\a b -> max a b) 0 $! R.map (\coeff -> normSquared coeff) array
            normSquared cn = (magnitude cn) P.^ 2
            norm :: Double
            norm = sqrt $! maxNormSquared

{-
-- Will have to encode positionRange in Position to be able to traverse index map
instance NormalizeVector (InfKet (Complex Double) Position) where
    normalizeVector ket@(InfKet coeffMap) = (1 / (norm :+ 0)) *^ ket
        where
            -- ketNormSquared = foldl' (P.+) 0 $! map (normSquared .coeffMap) positionRange
            -- normalizing with max of norms instead of sum of squares of norm to make the states more visible in simulation
            maxNormSquared = foldl' (\a b -> max a b) 0 $! map (normSquared .coeffMap) positionRange
            normSquared cn = (magnitude cn) P.^ 2
            norm :: Double
            norm = sqrt $! maxNormSquared
-}


  
-- Hermitian Inner Product
instance (InnerProductSpace (Complex Double) (SparseKet (Complex Double) Position)) where
    (<.>) (SparseKet coeffA) (SparseKet coeffB) = M.foldl' (+) 0 $! M.intersectionWith (\ampA ampB -> ampA * (complexConjugate ampB)) coeffA coeffB



-------------------
-- For field instance of Complex Double

instance Euclidean (Complex Double) where degree _ = Just 1
instance Division (Complex Double) where (/) = (P./)
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
instance Group (Complex Double) where (-) = (P.-)
instance LeftModule Natural (Complex Double) where
    n .* m = (fromIntegral $ toInteger n) P.* m
instance RightModule Natural (Complex Double) where
    m *. n = m P.* (fromIntegral n)
instance LeftModule Integer (Complex Double) where
    n .* m = (fromIntegral $ toInteger n) P.* m
instance RightModule Integer (Complex Double) where
    m *. n = m P.* (fromIntegral n)
instance Abelian (Complex Double)




--------------
-- For field instance of Double

instance Euclidean Double where degree _ = Just 1
instance Division Double where (/) = (P./)
instance PID Double
instance UFD Double
instance GCDDomain Double
instance IntegralDomain Double
instance Unital Double where one = 1
instance DecidableZero Double where isZero = (0 ==)
instance ZeroProductSemiring Double
instance Commutative Double
instance Multiplicative Double where (*) = (P.*)
instance DecidableUnits Double where recipUnit a = case a of
                                                       1 -> Just a
                                                       -1 -> Just a
                                                       _ -> Nothing
instance UnitNormalForm Double
instance Ring Double
instance Monoidal Double where zero = 0
instance Semiring Double
instance Additive Double where (+) = (P.+)
instance DecidableAssociates Double where isAssociate = (==)
instance Rig Double
instance Group Double where (-) = (P.-)
instance LeftModule Natural Double where
    n .* m = (fromIntegral $ toInteger n) P.* m
instance RightModule Natural Double where
    m *. n = m P.* (fromIntegral n)
instance LeftModule Integer Double where
    n .* m = (fromIntegral $ toInteger n) P.* m
instance RightModule Integer Double where
    m *. n = m P.* (fromIntegral n)
instance Abelian Double



---------------
-- For use by repa
instance Elt (Complex Double) where
    touch (a :+ b) = touch a >> touch b 
