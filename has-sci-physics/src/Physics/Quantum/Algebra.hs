{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies, FlexibleInstances, UndecidableInstances, GADTs, TypeOperators, FlexibleContexts #-}

module Physics.Quantum.Algebra where

import Numeric.Field.Class          (Field)
import Numeric.Algebra.Class        (LeftModule((.*)), Multiplicative((*)), Monoidal(zero))
import Numeric.Additive.Class       (Additive((+)))
import Numeric.Algebra.Unital       (Unital(one))

import Prelude hiding ((+), (*), Complex)

import qualified Data.Map as M
import qualified Data.List as L
import Data.Maybe (catMaybes)

import Data.Complex                 (Complex, conjugate)

import qualified Data.Array.Repa as R
import qualified Data.Array.Repa.Eval as R
import qualified Data.Array.Repa.Repr.Unboxed as R
import qualified Data.Array.Repa.Repr.Vector as R

import Data.Functor.Identity (runIdentity)

import qualified Data.Eigen.SparseMatrix as ES
import qualified Data.Eigen.Matrix as E





-- a vector space m over field k
class (Field k, Additive m) => VectorSpace k m | m -> k where
    (*^) :: k -> m -> m

{-
-- lhs type ‘m’ does not determine rhs type ‘k’
instance (Field k, LeftModule k m) => VectorSpace k m where
    (*^) = (.*)
-}

class VectorSpace k m => LieAlgebra k m | m -> k where
    lieBr :: m -> m -> m


class (VectorSpace k a, VectorSpace k b) => DualVectors k a b | a b -> k where
    dual :: a -> b

class VectorSpace k m => InnerProductSpace k m | m -> k where
    (<.>) :: m -> m -> k



-- Endo-morphism of vector space
-- Data structure for operator may be different for different use cases
-- A Map be suitable for sparse operator, matrix for discreet vector indices, and functions for continuous vector indices.
class VectorSpace k v => SimpleOperator k v o | o -> k, o -> v where
    op :: o -> v -> v

-- Once we have simple operator, we can extract the map information from that
-- and get addition and multiplication operations
data OperatorMap v where
    OperatorMap :: VectorSpace k v => (v -> v) -> OperatorMap v

unOperatorMap :: OperatorMap v -> (v -> v)
unOperatorMap (OperatorMap m) = m

extractOperatorMap :: SimpleOperator k v o => o -> OperatorMap v
extractOperatorMap = OperatorMap . op

instance Additive v => Additive (OperatorMap v) where
    (+) (OperatorMap operator1) (OperatorMap operator2) = OperatorMap $ \vec -> (operator1 vec) + (operator2 vec)

instance  Multiplicative (OperatorMap o) where
    (*) (OperatorMap operator1) (OperatorMap operator2) = OperatorMap $ \vec -> operator1 (operator2 vec)



class (DualVectors k cv v, SimpleOperator k v o) => Operator k v cv o | o -> k where
    (|><|) :: v -> cv -> o


(<|>) :: (InnerProductSpace k v, DualVectors k cv v) => cv -> v -> k
(<|>) covec vec = (dual covec) <.> vec


(|>) :: SimpleOperator k v o => o -> v -> v
(|>) = op


(<|) :: (InnerProductSpace k v, DualVectors k cv v, Operator k v cv o) =>  cv -> o -> (v -> k)
(<|) covector operator  = \vector -> covector <|> (operator |> vector)

    

class IsComplex a where
    complexConjugate :: a -> a

instance Num a => IsComplex (Complex a) where
    complexConjugate = conjugate

    
----------------------------------
-- Constructions of vector spaces


-- m is a set of vector indices indices
-- May or may not represent orthogonal vectors, write instance of InnerProductSpace separately
-- kets are represented as maps from indices to coefficients
-- sparse kets have finite number of non zero entries
data SparseKet k m where
    SparseKet :: (Field k, Ord m) => M.Map m k -> SparseKet k m

    
unSparseKet :: SparseKet k m -> M.Map m k
unSparseKet (SparseKet ketMap) = ketMap
    
data SparseBra k m where
     SparseBra :: (Field k, Ord m) => M.Map m k -> SparseBra k m

unSparseBra :: SparseBra k m -> M.Map m k
unSparseBra (SparseBra braMap) = braMap


instance (Show k, Show m) => Show (SparseKet k m) where
    show ket@(SparseKet map)  = show map

instance (Show k, Show m) => Show (SparseBra k m) where
    show ket@(SparseBra map)  = show map

instance (Ord m, Field k) => Additive (SparseKet k m) where
    (+) (SparseKet a) (SparseKet b) = SparseKet $ M.unionWith (+) a b

instance (Ord m, Field k) => Additive (SparseBra k m) where
    (+) (SparseBra a) (SparseBra b) = SparseBra $ M.unionWith (+) a b

instance (Ord m, Field k) => VectorSpace k (SparseKet k m) where
    (*^) s (SparseKet v) = SparseKet $ M.map (s *) v

instance (Field k, Ord m) => VectorSpace k (SparseBra k m) where
    (*^) s (SparseBra v) = SparseBra $ M.map (s *) v

instance (Field k, IsComplex k, Ord m) => DualVectors k (SparseKet k m) (SparseBra k m) where
    dual (SparseKet a) = SparseBra $ M.map complexConjugate a

instance (Field k, IsComplex k, Ord m) => DualVectors k (SparseBra k m) (SparseKet k m) where
    dual (SparseBra a) = SparseKet $ M.map complexConjugate a



-- very generic mappings
data InfKet k m = InfKet { unContKet :: (m -> k) }
data InfBra k m = InfBra { unContBra :: (m -> k) }


instance Field k => Additive (InfKet k m) where
    (+) (InfKet a) (InfKet b) = InfKet $ \x -> a x + b x

instance Field k => Additive (InfBra k m) where
    (+) (InfBra a) (InfBra b) = InfBra $ \x -> a x + b x

instance Field k => VectorSpace k (InfKet k m) where
    (*^) s (InfKet v) = InfKet $ \x -> s * (v x)
   
instance Field k => VectorSpace k (InfBra k m) where
    (*^) s (InfBra v) = InfBra $ \x -> s * (v x)

instance (Field k, IsComplex k) => DualVectors k (InfKet k m) (InfBra k m) where
    dual (InfKet a) = InfBra $ \x -> complexConjugate $ a x

instance (Field k, IsComplex k) => DualVectors k (InfBra k m) (InfKet k m) where
    dual (InfBra a) = InfKet $ \x -> complexConjugate $ a x



-- a vector of coefficients at lattice points
data LatticeKet k m = LatticeKet { unLatticeKet :: R.Array R.D R.DIM1 k }
data LatticeBra k m = LatticeBra { unLatticeBra :: R.Array R.D R.DIM1 k }


instance (Show k, Show m) => Show (LatticeKet k m) where
    -- Delayed representation doesn't have show instance
    show (LatticeKet vec) = show $! compute vec
        where compute :: Show k => R.Array R.D R.DIM1 k -> R.Array R.V R.DIM1 k
              compute = R.computeS

instance (Show k, Show m) => Show (LatticeBra k m) where
    -- Delayed representation doesn't have show instance
    show (LatticeBra vec) = show $! compute vec
        where compute :: Show k => R.Array R.D R.DIM1 k -> R.Array R.V R.DIM1 k
              compute = R.computeS

instance Field k => Additive (LatticeKet k m) where
   (+) (LatticeKet a) (LatticeKet b) = LatticeKet $! R.zipWith (+) a b

instance Field k => Additive (LatticeBra k m) where
   (+) (LatticeBra a) (LatticeBra b) = LatticeBra $! R.zipWith (+) a b

instance Field k => VectorSpace k (LatticeKet k m) where
    (*^) s (LatticeKet v) = LatticeKet $! R.map (s *) v
 
instance Field k => VectorSpace k (LatticeBra k m) where
    (*^) s (LatticeBra v) = LatticeBra $! R.map (s *) v

instance (Field k, IsComplex k) => DualVectors k (LatticeKet k m) (LatticeBra k m) where
    dual (LatticeKet a) = LatticeBra $! R.map complexConjugate a

instance (Field k, IsComplex k) => DualVectors k (LatticeBra k m) (LatticeKet k m) where
    dual (LatticeBra a) = LatticeKet $! R.map complexConjugate a



-- n x 1 matrix
data SparseMatKet k ck m where
    SparseMatKet :: E.Elem k ck => ES.SparseMatrix k ck -> SparseMatKet k ck m
    
unSparseMatKet :: SparseMatKet k ck m -> ES.SparseMatrix k ck
unSparseMatKet (SparseMatKet mat) = mat

data SparseMatBra k ck m where
    SparseMatBra :: E.Elem k ck => ES.SparseMatrix k ck -> SparseMatBra k ck m
    
unSparseMatBra :: SparseMatBra k ck m -> ES.SparseMatrix k ck
unSparseMatBra (SparseMatBra mat) = mat


instance (Show k) => Show (SparseMatKet k ck m) where
    show (SparseMatKet mat) = show mat

instance (Show k) => Show (SparseMatBra k ck m) where
    show (SparseMatBra mat) = show mat  

instance Additive (SparseMatKet k ck m) where
    (+) (SparseMatKet a) (SparseMatKet b) = SparseMatKet $! ES.add a b

instance Additive (SparseMatBra k ck m) where
    (+) (SparseMatBra a) (SparseMatBra b) = SparseMatBra $! ES.add a b

instance Field k => VectorSpace k (SparseMatKet k ck m) where
    (*^) s (SparseMatKet v) = SparseMatKet $! ES.scale s v
 
instance Field k => VectorSpace k (SparseMatBra k ck m) where
    (*^) s (SparseMatBra v) = SparseMatBra $! ES.scale s v


-------------------------------------
-- Constructions of operators


-- Operator formed from outerproduct of two vectors
data OuterOperator k v cv where
    OuterOperator :: (InnerProductSpace k v, DualVectors k cv v) => [(v, cv)] -> OuterOperator k v cv
    
unOuterOperator :: OuterOperator k v cv -> [(v, cv)]
unOuterOperator (OuterOperator outerProducts) = outerProducts

instance ( InnerProductSpace k v
         , DualVectors k cv v
         )
         => SimpleOperator k v (OuterOperator k v cv) where
    op (OuterOperator outerProducts) vec = L.foldr1 (+) $ map (\(ovec,ocovec) -> (ocovec <|> vec) *^ ovec) $ outerProducts

instance ( InnerProductSpace k v
         , DualVectors k cv v
         )
         => Operator k v cv (OuterOperator k v cv) where
    (|><|) vec covec = OuterOperator [(vec, covec)]

instance Additive (OuterOperator k v cv) where
    (+) (OuterOperator o1) (OuterOperator o2) = OuterOperator $ o1 L.++ o2




-- a sparse operator represented as a map from (i,j) to coefficient of |i><j|
-- Note that (|i><j|) |k> = <j|k> |i> and therefore the operation requires inner product to be defined
data SparseOperator k m = SparseOperator { unSparseOperator :: M.Map (m, m) k }


instance ( Field k, Ord m,
           InnerProductSpace k (SparseKet k m),
           DualVectors k (SparseBra k m) (SparseKet k m)
         )
         => SimpleOperator k (SparseKet k m) (SparseOperator k m) where  
    op (SparseOperator opMap) ket = SparseKet $! M.fromListWith (+) $!
        map (\(veci, coveci) -> ( veci
                                , (SparseBra $! M.fromListWith (+) $! catMaybes $! [fmap ((,) coveci) $! M.lookup (veci, coveci) opMap]) <|> ket
                                )
            ) $!
            (M.keys opMap) 

instance ( Field k, Ord m,
           InnerProductSpace k (SparseKet k m),
           DualVectors k (SparseBra k m) (SparseKet k m)
         )
         => Operator k (SparseKet k m) (SparseBra k m) (SparseOperator k m) where  
    (|><|) (SparseKet vecMap) (SparseBra covecMap) = SparseOperator $! M.fromListWith (+) $!
        catMaybes [(,) (i,j) <$> (fmap (*) (M.lookup i vecMap) <*> (M.lookup j covecMap)) | i <- M.keys vecMap, j <- M.keys covecMap ]

instance ( Field k, Ord m,
           InnerProductSpace k (SparseKet k m),
           DualVectors k (SparseBra k m) (SparseKet k m)
         )
         => Additive (SparseOperator k m) where  
    (+) (SparseOperator a) (SparseOperator b) = SparseOperator $! M.unionWith (+) a b

instance ( Field k, Ord m,
           InnerProductSpace k (SparseKet k m),
           DualVectors k (SparseBra k m) (SparseKet k m)
         )
         => Multiplicative (SparseOperator k m) where  
     (*) (SparseOperator a) (SparseOperator b) = SparseOperator $! M.fromListWith (+) $! multiplyMaps a b []
         where
             -- chaoose one index i and sum (over j and k) |j><i| |i><k|
             multiplyMaps :: ( Field k, Ord m,
                               InnerProductSpace k (SparseKet k m),
                               DualVectors k (SparseBra k m) (SparseKet k m)
                             )
                             => M.Map (m,m) k -> M.Map (m,m) k -> [((m,m), k)] -> [((m,m), k)]
             multiplyMaps mapA mapB oldResult = case null mapA || null mapB of
                 True -> oldResult
                 False ->
                   -- toList uses list fusion so this should not go through entire map
                   let oneElement = L.take 1 $ M.toList $ mapA
                   in case oneElement of
                          [] -> oldResult
                          elementToUse : _ -> let indexToUse = (snd $! fst elementToUse)

                                                  -- <i|i>
                                                  indexKet = SparseKet $! M.singleton indexToUse one
                                                  indexSquaredNorm =  indexKet <.> indexKet

                                                  -- take all keys from mapA which have this index in snd of key
                                                  (mapAToUse, mapAToPass) = M.partitionWithKey (\key _ -> snd key == indexToUse) mapA
                                                  -- take all keys from mapA which have this index in snd of key
                                                  (mapBToUse, mapBToPass) = M.partitionWithKey (\key _ -> fst key == indexToUse) mapB

                                                  -- we can discard second index from mapAToUse and first index from mapBToUse
                                                  -- we also know that the preserved indices will occur uniquely in index
                                                  -- so we can just take the cross product
                                                  listAToUse = M.toList $! mapAToUse
                                                  listBToUse = M.toList $! mapBToUse
                                                  newResult = [((a,b), va * vb * indexSquaredNorm) | ((a,_), va) <- listAToUse,
                                                                                                     ((_,b), vb) <- listBToUse
                                                                                                   ]

                                              in multiplyMaps mapAToPass mapBToPass $! newResult ++ oldResult

instance ( Field k, Ord m,
           InnerProductSpace k (SparseKet k m),
           DualVectors k (SparseBra k m) (SparseKet k m)
         )
         => LeftModule k (SparseOperator k m) where
    (.*) scalar (SparseOperator b) = SparseOperator $! M.map (scalar *) b



-- a more general operator represented as a map from j to Map i_n s_n where operator takes |i_n> to s_n |j>
-- this construct may not suffice when there are uncountable such i_n
data InfOperator k m = InfOperator { unInfOperator :: m -> M.Map m k }

instance (Field k) => SimpleOperator k (InfKet k m) (InfOperator k m) where
    op operator = \(InfKet coeffMap) -> InfKet $ \index -> foldr1 (+) $
        map ( \(invIndex, factor) -> factor * coeffMap invIndex ) $
            M.toList $ unInfOperator operator $ index

instance (Field k, Ord m) => Additive (InfOperator k m) where
    (+) (InfOperator op1) (InfOperator op2) = InfOperator $ \index -> M.unionWith (+) (op1 index) (op2 index)

instance (Field k, Ord m) => Multiplicative (InfOperator k m) where
    (*) (InfOperator op1) (InfOperator op2) = InfOperator $ \index ->
            foldr1 (M.unionWith (+)) $
                map (\(i,s) -> M.map (s *) $ op2 i) $
                    M.toList $ op1 index



-- A matrix to transfom LatticeKet vector
data LatticeOperator k m where
    LatticeOperator :: (Field k, R.Elt k, R.Unbox k, Num k) => R.Array R.D R.DIM2 k -> LatticeOperator k m

unLatticeOperator :: LatticeOperator k m -> R.Array R.D R.DIM2 k
unLatticeOperator (LatticeOperator mat) = mat

-- Num k is required because of sumP in method definition
instance Field k => SimpleOperator k (LatticeKet k m) (LatticeOperator k m) where
    op (LatticeOperator matrix) (LatticeKet vector) = LatticeKet $! R.delay $! runIdentity $! mvMultP matrix vector

instance Additive (LatticeOperator k m) where
    (+) (LatticeOperator op1) (LatticeOperator op2) = LatticeOperator $! R.zipWith (+) op1 op2

-- Num k is required because of sumP in method definition
instance Multiplicative (LatticeOperator k m) where
    (*) (LatticeOperator op1) (LatticeOperator op2) = LatticeOperator $! R.delay $! runIdentity $! mmMultP op1 op2

instance Field k => LeftModule k (LatticeOperator k m) where
    (.*) scalar (LatticeOperator matrix) = LatticeOperator $! R.map (scalar *) matrix


-- Taken from Repa wiki, extended it to Additive and Multiplicative structures
mmMultP :: (Monad m, R.Source r k, Additive k, Multiplicative k, R.Elt k, R.Unbox k, Num k)
        => R.Array r R.DIM2 k
        -> R.Array r R.DIM2 k
        -> m (R.Array R.U R.DIM2 k)
mmMultP a b = R.sumP (R.zipWith (*) aRepl bRepl)
    where
      -- Note: transpose transposes the lowest two dimensions
      t = R.transpose b
      aRepl = R.extend (R.Z R.:.R.All R.:.colsB R.:.R.All) a
      bRepl = R.extend (R.Z R.:.rowsA R.:.R.All R.:.R.All) t

      (R.Z R.:.colsA R.:.rowsA) = R.extent a
      (R.Z R.:.colsB R.:.rowsB) = R.extent b

-- Multiply a martix with a vector
-- based on the above matrix matrix multiplication
mvMultP :: (Monad m, R.Source r k, Additive k, Multiplicative k, R.Elt k, R.Unbox k, Num k)
        => R.Array r R.DIM2 k
        -> R.Array r R.DIM1 k
        -> m (R.Array R.U R.DIM1 k)
mvMultP a b = R.sumP (R.zipWith (*) aRepl bRepl)
    where
      aRepl = a
      bRepl = R.extend (R.Z R.:.rowsA R.:.R.All) b
      (R.Z R.:.rowsA R.:.colsA) = R.extent a




-- A matrix to transfom LatticeKet vector
data SparseMatOperator k ck m where
    SparseMatOperator :: E.Elem k ck => ES.SparseMatrix k ck -> SparseMatOperator k ck m
    
unSparseMatOperator :: SparseMatOperator k ck m -> ES.SparseMatrix k ck
unSparseMatOperator (SparseMatOperator mat) = mat


-- Num k is required because of sumP in method definition
instance Field k => SimpleOperator k (SparseMatKet k ck m) (SparseMatOperator k ck m) where
    op (SparseMatOperator matrix) (SparseMatKet vector) = SparseMatKet $! ES.mul matrix vector

instance Field k => Additive (SparseMatOperator k ck m) where
    (+) (SparseMatOperator op1) (SparseMatOperator op2) = SparseMatOperator $! ES.add op1 op2

instance Field k => Multiplicative (SparseMatOperator k ck m) where
    (*) (SparseMatOperator op1) (SparseMatOperator op2) = SparseMatOperator $! ES.mul op1 op2

instance Field k => LeftModule k (SparseMatOperator k ck m) where
    (.*) scalar (SparseMatOperator matrix) = SparseMatOperator $! ES.scale scalar matrix
