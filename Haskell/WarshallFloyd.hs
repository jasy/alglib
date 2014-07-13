#!/usr/bin/env runghc
import Control.Monad
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Unboxed.Mutable as VM

wf :: Int -> [(Int,Int,Int)] -> V.Vector (VU.Vector Int)
wf n es = V.generate n (\i->VU.slice (n*i) n d) where
  idx i j = i*n+j
  d = VU.create $ do
    d <- VM.replicate (n*n) (maxBound `div` 2)
    forM_ [0..n-1] $ \i -> VM.write d (idx i i) 0
    forM_ es $ \(i, j, t) -> do
      VM.write d (idx i j) t
      VM.write d (idx j i) t
    forM_ [0..n-1] $ \k ->
      forM_ [0..n-1] $ \i ->
        forM_ [0..n-1] $ \j -> do
          a <- VM.read d (idx i j)
          b <- liftM2 (+) (VM.read d (idx i k)) (VM.read d (idx k j))
          VM.write d (idx i j) $ min a b
    return d
