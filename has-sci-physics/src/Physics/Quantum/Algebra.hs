{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies, FlexibleInstances, UndecidableInstances, GADTs #-}

module Physics.Quantum.Algebra where

import Numeric.Field.Class          (Field)
import Numeric.Algebra.Class        (LeftModule((.*)), Multiplicative((*)))
import Numeric.Additive.Class       (Additive(..))


import Prelude hiding ((+), (*))

import qualified Data.Map as M
import Data.Maybe (catMaybes)



-- a vector space m over field k
class (Field k, Additive m) => VectorSpace k m | m -> k where
    (*^) :: k -> m -> m


{- The LeftModule condition breaks the functional dependency of vector space
instance (Additive m, Field k, LeftModule k m) => VectorSpace k m where
    (*^) = (.*)
-}

class VectorSpace k m => LieAlgebra k m | m -> k where
    lieBr :: m -> m -> m


class (VectorSpace k a, VectorSpace k b) => DualVectors k a b | a b -> k where
    dual :: a -> b

class VectorSpace k m => InnerProductSpace k m | m -> k where
    (<.>) :: m -> m -> k

class VectorSpace k m => OuterProductSpace k m | m -> k where
    (>.<) :: m -> m -> (m -> m -> k)

-- For some reason it doesn't deduct VectorSpace contraint from InnerProductSpace constraint
instance (VectorSpace k m, InnerProductSpace k m)  => OuterProductSpace k m where
    (>.<) a b = \c d ->  (a <.> c) * (b <.> d)


-- Endo-morphism of vector space
-- Data structure for operator may be different for different use cases
-- A Map be suitable for sparse operator, matrix for discreet vector indices, and functions for continuous vector indices.
class VectorSpace k v => SimpleOperator k v o | o -> k where
    op :: o -> v -> v

class (DualVectors k cv v, SimpleOperator k v o) => Operator k v cv o | o -> k where
    (|><|) :: v -> cv -> o
    

(<|>) :: (InnerProductSpace k v, DualVectors k cv v) => cv -> v -> k
(<|>) covec vec = (dual covec) <.> vec


(|>) :: SimpleOperator k v o => o -> v -> v
(|>) = op


(<|) :: (InnerProductSpace k v, DualVectors k cv v, Operator k v cv o) => o -> cv -> (v -> k)
(<|) operator covector = \vector -> covector <|> (operator |> vector)


class Complex a where
    complexConjugate :: a -> a


------ 
-- Constructions of vector spaces


-- m is a set of vector indices indices
-- May or may not represent orthogonal vectors, write instance of InnerProductSpace separately
-- kets are represented as maps from indices to coefficients
-- sparse kets have finite number of non zero entries
data SparseKet k m where
    SparseKet :: (Field k, Ord m) => M.Map m k -> SparseKet k m
    
unSparseKet :: SparseKet k m -> M.Map m k
unSparseKet (SparseKet ketMap) = ketMap

  
data InfKet k m = InfKet { unContKet :: (m -> k) }
    
data SparseBra k m where
     SparseBra :: (Field k, Ord m) => M.Map m k -> SparseBra k m

unSparseBra :: SparseBra k m -> M.Map m k
unSparseBra (SparseBra braMap) = braMap

   
data InfBra k m = InfBra { unContBra :: (m -> k) }

data SparseOperator k m = SparseOperator { unSparseOperator :: M.Map (m, m) k }
data InfOperator k m = InfOperator { unInfOperator :: m -> m -> k }


instance ( Field k, Ord m,
           InnerProductSpace k (SparseKet k m),
           DualVectors k (SparseBra k m) (SparseKet k m)
         )
         => SimpleOperator k (SparseKet k m) (SparseOperator k m) where
  
    op (SparseOperator opMap) ket = SparseKet $ M.fromListWith (+) $
        map (\(veci, coveci) -> (veci, (SparseBra $ M.fromList $ catMaybes $ [fmap ((,) coveci) $ M.lookup (veci, coveci) opMap] ) <|> ket ) ) $
            (M.keys opMap)
 


instance ( Field k, Ord m,
           InnerProductSpace k (SparseKet k m),
           DualVectors k (SparseBra k m) (SparseKet k m)
         )
         => Operator k (SparseKet k m) (SparseBra k m) (SparseOperator k m) where
  
    (|><|) (SparseKet vecMap) (SparseBra covecMap) = SparseOperator $ M.fromList $
        catMaybes [(,) (i,j) <$> (fmap (*) (M.lookup i vecMap) <*> (M.lookup j covecMap)) | i <- M.keys vecMap, j <- M.keys covecMap ]



instance (Ord m, Field k) => Additive (SparseKet k m) where
    (+) (SparseKet a) (SparseKet b) = SparseKet $ M.unionWith (+) a b

instance Field k => Additive (InfKet k m) where
    (+) (InfKet a) (InfKet b) = InfKet $ \x -> a x + b x

instance (Ord m, Field k) => Additive (SparseBra k m) where
    (+) (SparseBra a) (SparseBra b) = SparseBra $ M.unionWith (+) a b

instance Field k => Additive (InfBra k m) where
    (+) (InfBra a) (InfBra b) = InfBra $ \x -> a x + b x


    
instance (Ord m, Field k) => VectorSpace k (SparseKet k m) where
    (*^) s (SparseKet v) = SparseKet $ M.map (s *) v
    
instance Field k => VectorSpace k (InfKet k m) where
    (*^) s (InfKet v) = InfKet $ \x -> s * (v x)

instance (Field k, Ord m) => VectorSpace k (SparseBra k m) where
    (*^) s (SparseBra v) = SparseBra $ M.map (s *) v
    
instance Field k => VectorSpace k (InfBra k m) where
    (*^) s (InfBra v) = InfBra $ \x -> s * (v x)    



instance (Field k, Complex k, Ord m) => DualVectors k (SparseKet k m) (SparseBra k m) where
    dual (SparseKet a) = SparseBra $ M.map complexConjugate a

instance (Field k, Complex k) => DualVectors k (InfKet k m) (InfBra k m) where
    dual (InfKet a) = InfBra $ \x -> complexConjugate $ a x

instance (Field k, Complex k, Ord m) => DualVectors k (SparseBra k m) (SparseKet k m) where
    dual (SparseBra a) = SparseKet $ M.map complexConjugate a

instance (Field k, Complex k) => DualVectors k (InfBra k m) (InfKet k m) where
    dual (InfBra a) = InfKet $ \x -> complexConjugate $ a x
