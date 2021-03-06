{-# LANGUAGE TypeOperators, DataKinds, KindSignatures, MultiParamTypeClasses, GADTs, StandaloneDeriving #-}


module Physics.Classical where

import qualified Numeric.Units.Dimensional as D
import qualified Numeric.Units.Dimensional.Quantities as D
import Numeric.Units.Dimensional                    ((*~))

import Numeric.Units.Dimensional.SIUnits            (newton, meter, kilo, gram, second, coulomb)
import Numeric.NumType.DK.Integers                  (pos2, pos3, TypeInt(Pos2, Pos3))



-----------------
-- Rewrite this module for n-dimensional vectors, something like in Physics.Quantum.Algebra
-- Add electromagnetism and rotation
-----------------

data Vec (b::D.Dimension) a where
    Vec :: (Eq a, Floating a) => (D.Quantity b a, D.Quantity b a) -> Vec b a

deriving instance Eq (Vec b a)

getVec :: Vec b a -> (D.Quantity b a, D.Quantity b a)
getVec (Vec vec) = vec


-- vector addition
(^+) :: Vec b a -> Vec b a -> Vec b a
(^+) (Vec (x,y)) (Vec (l,m)) = Vec (x D.+ l, y D.+ m)

-- vector subtraction
(^-) :: Vec b a -> Vec b a -> Vec b a
(^-) (Vec (x,y)) (Vec (l,m)) = Vec (x D.- l, y D.- m)

-- scalar product with a product
(*^) :: D.Quantity b a -> Vec c a -> Vec (b D.* c) a
(*^) x (Vec (l,m)) = Vec (x D.* l, x D.* m)

-- scalar (dot) product of two vectors
(.<) :: Vec b a -> Vec c a -> D.Quantity (b D.* c) a
(.<) (Vec (x,y)) (Vec (l,m)) = (x D.* l D.+ y D.* m)


norm :: Floating a => Vec b a -> D.Quantity (D.Root (b D.* b) Pos2) a
norm vec = D.sqrt (vec .< vec)

-- zeroVector :: Num a => Vec b a
-- zeroVector = Vec (_0, _0)


type SpatialField a b = Vec D.DLength a -> b

type GravityFieldUnit = D.DLength D./ D.DTime D.^ Pos2
type ElectricFieldUnit = D.DForce D./ D.DElectricCharge

class (Eq b, Floating b) => HasMass a b where
    mass :: a -> D.Mass b

class (Eq b, Floating b) => HasElectricCharge a b where
    electricCharge :: a -> D.Quantity D.DElectricCharge b

class HasCenterOfGravity a b where
    centerOfGravity :: a -> Vec D.DLength b

class ModifiableCenterOfGravity a b where
    changeCenterOfGravity :: Vec D.DLength b -> a -> a

class HasVelocity a b where
    velocity :: a -> Vec (D.DLength D./ D.DTime) b

class ModifiableVelocity a b where
    changeVelocity :: Vec (D.DLength D./ D.DTime) b -> a -> a




gravitationalConstant :: Floating a => D.Quantity (D.DLength D.^ Pos3  D./ (D.DMass D.* D.DTime D.^ Pos2)) a
gravitationalConstant = 6.673e-11 *~ (newton D.* meter D.^ pos2 D./ (kilo gram) D.^ pos2)

electricFieldConstant :: Floating a => D.Quantity (D.DForce D.* D.DLength D.^ Pos2  D./ (D.DElectricCharge D.^ Pos2)) a
electricFieldConstant = 8.9875517e9 *~ (newton D.* meter D.^ pos2 D./ coulomb D.^ pos2)

nBodyGravitationField :: (HasMass a e, HasCenterOfGravity a e) => [a] -> SpatialField e (Vec GravityFieldUnit e)
nBodyGravitationField particles =
    let gravityByOneParticle = \particle -> (\x -> if x == centerOfGravity particle
                                                   then zeroGravityValue
                                                   else ( gravitationalConstant D.* (mass particle) D./ ( (norm $ centerOfGravity particle ^- x) D.^ pos3 ) ) *^
                                                        (centerOfGravity particle ^- x)
                                            )        
    in foldr (\x y -> (\z -> x z ^+ y z))
             (\x -> zeroGravityValue) $
             map gravityByOneParticle particles

    where
        zeroGravityValue = Vec (0 *~ (meter D./ second D.^ pos2) , 0 *~ (meter D./ second D.^ pos2))
             


nBodyGravityStep :: ( HasMass a b, HasCenterOfGravity a b, HasVelocity a b
                    , ModifiableVelocity a b, ModifiableCenterOfGravity a b
                    )
                    => [a] -> D.Time b -> [a]
nBodyGravityStep initStates timeStep =
    let combinedGravityField = nBodyGravitationField initStates
    in map (\state -> changeVelocity (velocity state ^+ ( timeStep *^ (combinedGravityField $ centerOfGravity state) )) $
                      changeCenterOfGravity (centerOfGravity state ^+ ( timeStep *^ (velocity state) )) $ state
           ) initStates



nBodyElectricField :: (HasElectricCharge a e, HasCenterOfGravity a e) => [a] -> SpatialField e (Vec ElectricFieldUnit e)
nBodyElectricField particles =
    let electricFieldByOneParticle = \particle -> (\x -> if x == centerOfGravity particle
                                                   then zeroElectricFieldValue
                                                   else ( electricFieldConstant D.* (electricCharge particle) D./ ( (norm $ centerOfGravity particle ^- x) D.^ pos3 ) ) *^
                                                        (x ^- centerOfGravity particle)
                                            )        
    in foldr (\x y -> (\z -> x z ^+ y z))
             (\x -> zeroElectricFieldValue) $
             map electricFieldByOneParticle particles

    where
        zeroElectricFieldValue = Vec (0 *~ (newton D./ coulomb) , 0 *~ (newton D./ coulomb))


nBodyElectricPotential :: (HasElectricCharge a e, HasCenterOfGravity a e) => [a] -> SpatialField e (D.Quantity (ElectricFieldUnit D.* D.DLength) e)
nBodyElectricPotential particles =
    let electricPotentialByOneParticle = \particle -> (\x -> if x == centerOfGravity particle
                                                             then D._0
                                                             else ( electricFieldConstant D.* (electricCharge particle) D./ ( norm $ centerOfGravity particle ^- x) )

                                                      )        
    in foldr (\x y -> (\z -> x z D.+ y z))
             (\x -> D._0) $
             map electricPotentialByOneParticle particles



nBodyElectricStep :: ( HasElectricCharge a b, HasMass a b, HasCenterOfGravity a b, HasVelocity a b
                     , ModifiableVelocity a b, ModifiableCenterOfGravity a b
                     )
                    => [a] -> D.Time b -> [a]
nBodyElectricStep initStates timeStep =
    let combinedElectricField = nBodyElectricField initStates
    in map (\state -> changeVelocity (velocity state ^+ ( timeStep *^ ( ((electricCharge state) D./ (mass state)) *^
                                                                        (combinedElectricField $ centerOfGravity state)
                                                                      )
                                                        )
                                     ) $
                                     changeCenterOfGravity (centerOfGravity state ^+ ( timeStep *^ (velocity state) )) $ state
           ) initStates
